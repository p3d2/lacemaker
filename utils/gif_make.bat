@echo off

:: Variables
set "input_folder=Z:\Desktop\render"
set "output_file=Z:\Desktop\3024_1to4.gif"

:: Clear output file if it already exists
if exist "%output_file%" (
  del "%output_file%"
)

:: Temporary folder for processed images
set "temp_folder=%input_folder%\temp"
mkdir "%temp_folder%"

for %%I in ("%input_folder%\*.png") do (

	rem Step 1: Make non-transparent pixels black
    magick "%%I" -colorspace sRGB -fill black +opaque none "%temp_folder%\%%~nI_black.png"
    
    rem Step 2: Expand black pixels by 2px, ensuring color space is set to sRGB to match
    magick "%temp_folder%\%%~nI_black.png" -colorspace sRGB -morphology Dilate Disk:2 "%temp_folder%\%%~nI_expanded.png"
    
    rem Step 3: Ensure both images are in the same color space before compositing
    magick "%temp_folder%\%%~nI_expanded.png" -colorspace sRGB "%%I" -colorspace sRGB -compose Over -composite "%temp_folder%\%%~nI_final.png"
	
	del "%temp_folder%\%%~nI_black.png"
	del "%temp_folder%\%%~nI_expanded.png"
)

:: Get the last image file name
dir /b /o:n /a:-d "%temp_folder%\*.png" > list.txt
for /f %%I in (list.txt) do set "last=%%I"
del list.txt

:: Convert images to GIF
setlocal enabledelayedexpansion
set "file_list="
for %%I in ("%temp_folder%\*.png") do (
  set "file_list=!file_list! "%%I""
)
magick convert -dispose previous -loop 0 -delay 10 !file_list! -delay 100 "%temp_folder%\%last%" "%output_file%"
endlocal

:: Clean up temporary folder
rmdir /s /q "%temp_folder%"
