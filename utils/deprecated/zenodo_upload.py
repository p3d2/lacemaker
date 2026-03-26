import os
import zipfile
import requests
from pathlib import Path
from tqdm.notebook import tqdm

# ── Configuration ─────────────────────────────────────────────────────────────
FOLDER     = Path(r"C:\path\to\your\folder")  # folder to zip and upload
ZIP_NAME   = "textile_samples"                # base name for output parts
RECORD_ID  = "19232529"
TOKEN      = "your_token_here"
CHUNK_GB   = 2.0    # size of each part in GB (keep under 4 GB for FAT32)
KEEP_PARTS = False  # set to True to keep part files after upload
ZENODO_API = "https://zenodo.org/api"
# ──────────────────────────────────────────────────────────────────────────────


def get_existing_files(record_id, token):
    url = f"{ZENODO_API}/records/{record_id}/draft/files"
    r = requests.get(url, headers={"Authorization": f"Bearer {token}"})
    r.raise_for_status()
    return {entry["key"] for entry in r.json()["entries"]}


def zip_folder(folder_path, zip_path):
    folder_path = Path(folder_path)
    files = sorted(f for f in folder_path.rglob("*") if f.is_file())
    total = sum(f.stat().st_size for f in files)
    print(f"  Zipping {len(files)} files ({total / 1e9:.2f} GB uncompressed)...")
    with tqdm(total=total, unit="B", unit_scale=True) as pbar:
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for file in files:
                zf.write(file, file.relative_to(folder_path.parent))
                pbar.update(file.stat().st_size)


def split_zip(zip_path, base_name, chunk_bytes):
    zip_size = zip_path.stat().st_size
    n_parts = (zip_size + chunk_bytes - 1) // chunk_bytes
    print(f"  Splitting into {n_parts} part(s) of {chunk_bytes / 1e9:.1f} GB...")
    parts = []
    with open(zip_path, "rb") as f:
        for i in range(n_parts):
            part_path = zip_path.parent / f"{base_name}.part{i+1:03d}"
            with open(part_path, "wb") as out:
                out.write(f.read(chunk_bytes))
            parts.append(part_path)
            print(f"    {part_path.name}  ({part_path.stat().st_size / 1e9:.2f} GB)")
    return parts


def upload_file(file_path, record_id, token):
    filename = file_path.name
    file_size = file_path.stat().st_size

    # 1. Initialise file slot
    r = requests.post(
        f"{ZENODO_API}/records/{record_id}/draft/files",
        headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"},
        json=[{"key": filename}],
    )
    r.raise_for_status()

    # 2. Stream upload
    with open(file_path, "rb") as f:
        with tqdm.wrapattr(f, "read", total=file_size, unit="B",
                           unit_scale=True, desc=f"  Uploading {filename}") as wrapped:
            r = requests.put(
                f"{ZENODO_API}/records/{record_id}/draft/files/{filename}/content",
                headers={"Authorization": f"Bearer {token}",
                         "Content-Type": "application/octet-stream"},
                data=wrapped,
            )
    r.raise_for_status()

    # 3. Commit
    r = requests.post(
        f"{ZENODO_API}/records/{record_id}/draft/files/{filename}/commit",
        headers={"Authorization": f"Bearer {token}"},
    )
    r.raise_for_status()


# ── Run ───────────────────────────────────────────────────────────────────────
zip_path = FOLDER.parent / f"{ZIP_NAME}.zip"
chunk_bytes = int(CHUNK_GB * 1e9)

print(f"Checking existing files on record {RECORD_ID}...")
existing = get_existing_files(RECORD_ID, TOKEN)
if existing:
    print(f"  Already uploaded: {', '.join(sorted(existing))}")

# Step 1: zip
if zip_path.exists():
    print(f"  {zip_path.name} already exists, skipping zip step.")
else:
    zip_folder(FOLDER, zip_path)
    print(f"  Zip size: {zip_path.stat().st_size / 1e9:.2f} GB")

# Step 2: split
parts = sorted(zip_path.parent.glob(f"{ZIP_NAME}.part*"))
if parts:
    print(f"  Found {len(parts)} existing part(s), skipping split step.")
else:
    parts = split_zip(zip_path, ZIP_NAME, chunk_bytes)

# Step 3: upload parts not yet on Zenodo
to_upload = [p for p in parts if p.name not in existing]
print(f"\nUploading {len(to_upload)} part(s)...")
for i, part in enumerate(to_upload, 1):
    print(f"\n[{i}/{len(to_upload)}] {part.name}")
    try:
        upload_file(part, RECORD_ID, TOKEN)
        print(f"  Done.")
    except requests.HTTPError as e:
        print(f"  ERROR: {e.response.status_code} — {e.response.text}")
        raise
    finally:
        if not KEEP_PARTS:
            part.unlink()

if not KEEP_PARTS:
    zip_path.unlink(missing_ok=True)

print(f"\nAll done. Publish at: https://zenodo.org/uploads/{RECORD_ID}")
print(f"\nTo reassemble and extract on Linux/Mac:")
print(f"  cat {ZIP_NAME}.part* > {ZIP_NAME}.zip && unzip {ZIP_NAME}.zip")
print(f"To reassemble on Windows:")
print(f"  copy /b {ZIP_NAME}.part* {ZIP_NAME}.zip")
