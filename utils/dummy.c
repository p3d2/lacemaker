#include <stdio.h>
#include <string.h>

void eglGetDisplay() {
    printf("Dummy eglGetDisplay called.\n");
}

void eglInitialize() {
    printf("Dummy eglInitialize called.\n");
}

void eglTerminate() {
    printf("Dummy eglTerminate called.\n");
}

void eglChooseConfig() {
    printf("Dummy eglChooseConfig called.\n");
}

void eglCreateContext() {
    printf("Dummy eglCreateContext called.\n");
}

void eglDestroyContext() {
    printf("Dummy eglDestroyContext called.\n");
}

void eglCreateWindowSurface() {
    printf("Dummy eglCreateWindowSurface called.\n");
}

void eglDestroySurface() {
    printf("Dummy eglDestroySurface called.\n");
}

void eglMakeCurrent() {
    printf("Dummy eglMakeCurrent called.\n");
}

void eglSwapBuffers() {
    printf("Dummy eglSwapBuffers called.\n");
}

void eglGetProcAddress() {
    printf("Dummy eglGetProcAddress called.\n");
}

const char* eglQueryString(void* dpy, int name) {
    printf("Dummy eglQueryString called.\n");
    return "Mock EGL String";
}

void* eglGetCurrentDisplay() {
    printf("Dummy eglGetCurrentDisplay called.\n");
    return NULL;
}

void* eglGetCurrentSurface() {
    printf("Dummy eglGetCurrentSurface called.\n");
    return NULL;
}

void* eglQueryContext() {
    printf("Dummy eglQueryContext called.\n");
    return NULL;
}

void* eglBindAPI() {
    printf("Dummy eglBindAPI called.\n");
    return NULL;
}

void* eglGetCurrentContext() {
    printf("Dummy eglGetCurrentContext called.\n");
    return NULL;
}

void* eglGetConfigs() {
    printf("Dummy eglGetConfigs called.\n");
    return NULL;
}

void* eglCreatePbufferSurface() {
    printf("Dummy eglCreatePbufferSurface called.\n");
    return NULL;
}

int eglGetConfigAttrib(void* dpy, void* config, int attribute, int* value) {
    printf("Dummy eglGetConfigAttrib called.\n");
    *value = 0;
    return 1; // Return EGL_TRUE
}

int eglGetError(void* dpy, void* config, int attribute, int* value) {
    printf("Dummy eglGetError called.\n");
    *value = 0;
    return 1; // Return EGL_TRUE
}

int eglSwapInterval(void* dpy, void* config, int attribute, int* value) {
    printf("Dummy eglSwapInterval called.\n");
    *value = 0;
    return 1; // Return EGL_TRUE
}

void* glXChooseFBConfig(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glXChooseFBConfig called.\n");
    return NULL;
}

void* glXGetConfig(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glXGetConfig called.\n");
    return NULL;
}

void* glXGetFBConfigAttrib(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glXGetFBConfigAttrib called.\n");
    return NULL;
}

void* glXGetVisualFromFBConfig(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glXGetVisualFromFBConfig called.\n");
    return NULL;
}

void* glGetString(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glGetString called.\n");
    return NULL;
}

void* glGetIntegerv(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glGetIntegerv called.\n");
    return NULL;
}

void* glXChooseVisual(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glXChooseVisual called.\n");
    return NULL;
}

void* glMatrixMode(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glMatrixMode called.\n");
    return NULL;
}

void* glLoadIdentity(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glLoadIdentity called.\n");
    return NULL;
}

void* glOrtho(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glOrtho called.\n");
    return NULL;
}

void* glGenTextures(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glGenTextures called.\n");
    return NULL;
}

void* glViewport(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glViewport called.\n");
    return NULL;
}

void* glTexParameterf(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glTexParameterf called.\n");
    return NULL;
}

void* glEnable(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glEnable called.\n");
    return NULL;
}

void* glClear(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glClear called.\n");
    return NULL;
}

void* glDisable(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glDisable called.\n");
    return NULL;
}

void* glTexSubImage2D(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glTexSubImage2D called.\n");
    return NULL;
}

void* glPixelStorei(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glPixelStorei called.\n");
    return NULL;
}

void* glDeleteTextures(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glDeleteTextures called.\n");
    return NULL;
}

void* glTexImage2D(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glTexImage2D called.\n");
    return NULL;
}

void* glBlendFunc(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glBlendFunc called.\n");
    return NULL;
}

void* glLoadMatrixf(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glLoadMatrixf called.\n");
    return NULL;
}

void* glBindTexture(void* dpy, int screen, const int *attrib_list, int *nelements) {
    printf("Dummy glBindTexture called.\n");
    return NULL;
}