// Copyright 2011 Google Inc. All Rights Reserved.
//
// This code is licensed under the same terms as WebM:
//  Software License Agreement:  http://www.webmproject.org/license/software/
//  Additional IP Rights Grant:  http://www.webmproject.org/license/additional/
// -----------------------------------------------------------------------------
//
//  Simple WebP file viewer.
//
// Compiling on linux:
//   sudo apt-get install libglut3-dev mesa-common-dev
//   gcc -o vwebp vwebp.c -O3 -lwebp -lglut -lGL
// Compiling on Mac + XCode:
//   gcc -o vwebp vwebp.c -lwebp -framework GLUT -framework OpenGL
//
// Author: Skal (pascal.massimino@gmail.com)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "webp/decode.h"
#include "webp/mux.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#ifdef FREEGLUT
#include <GL/freeglut.h>
#endif
#endif

#ifdef _MSC_VER
#define snprintf _snprintf
#endif

// Unfortunate global variables
static const WebPDecBuffer* kPic = NULL;
static const char* file_name = NULL;
static int print_info = 0;
static int exiting = 0;

//------------------------------------------------------------------------------
// Callbacks

static void HandleKey(unsigned char key, int pos_x, int pos_y) {
  (void)pos_x;
  (void)pos_y;
  if (key == 'q' || key == 'Q' || key == 27 /* Esc */) {
#ifdef FREEGLUT
    glutLeaveMainLoop();
    exiting = 1;
#else
    WebPFreeDecBuffer((WebPDecBuffer*)kPic);
    kPic = NULL;
    exit(0);
#endif
  } else if (key == 'i') {
    print_info = 1 - print_info;
    glutPostRedisplay();
  }
}

static void HandleReshape(int width, int height) {
  // TODO(skal): proper handling of resize, esp. for large pictures.
  // + key control of the zoom.
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

static void PrintString(const char* const text) {
  void* const font = GLUT_BITMAP_9_BY_15;
  int i;
  for (i = 0; text[i]; ++i) {
    glutBitmapCharacter(font, text[i]);
  }
}

static void HandleDisplay(void) {
  if (kPic == NULL) return;
  glClear(GL_COLOR_BUFFER_BIT);
  glPushMatrix();
  glPixelZoom(1, -1);
  glRasterPos2f(-1, 1);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, kPic->u.RGBA.stride / 4);
  glDrawPixels(kPic->width, kPic->height, GL_RGBA, GL_UNSIGNED_BYTE,
               (GLvoid*)kPic->u.RGBA.rgba);
  if (print_info) {
    char tmp[32];

    glColor4f(0.0, 0.0, 0.0, 0.0);
    glRasterPos2f(-0.95f, 0.90f);
    PrintString(file_name);

    snprintf(tmp, sizeof(tmp), "Dimension:%d x %d", kPic->width, kPic->height);
    glColor4f(0.0, 0.0, 0.0, 0.0);
    glRasterPos2f(-0.95f, 0.80f);
    PrintString(tmp);
  }
  glFlush();
}

static void Show(const WebPDecBuffer* const pic) {
  glutInitDisplayMode(GLUT_RGBA);
  glutInitWindowSize(pic->width, pic->height);
  glutCreateWindow("WebP viewer");
  glutReshapeFunc(HandleReshape);
  glutDisplayFunc(HandleDisplay);
  glutIdleFunc(NULL);
  glutKeyboardFunc(HandleKey);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  HandleReshape(pic->width, pic->height);
}

//------------------------------------------------------------------------------
// File decoding

static void* ReadFile(const char* const in_file, size_t* size) {
  int ok;
  size_t data_size = 0;
  void* data = NULL;
  FILE* const in = fopen(in_file, "rb");
  *size = 0;

  if (!in) {
    fprintf(stderr, "cannot open input file '%s'\n", in_file);
    return NULL;
  }
  fseek(in, 0, SEEK_END);
  data_size = ftell(in);
  fseek(in, 0, SEEK_SET);
  data = malloc(data_size);
  if (data == NULL) return 0;
  ok = (fread(data, data_size, 1, in) == 1);
  fclose(in);

  if (!ok) {
    fprintf(stderr, "Could not read %zu bytes of data from file %s\n",
            data_size, in_file);
    free(data);
    return NULL;
  }

  *size = data_size;
  return data;
}

static int Decode(WebPMux* const mux, const int frame_number,
                  WebPDecoderConfig* const config,
                  uint32_t* const duration) {
  WebPData image;
  uint32_t x_off = 0, y_off = 0;
  WebPDecBuffer* const output_buffer = &config->output;
  int ok = 0;

  WebPFreeDecBuffer((WebPDecBuffer*)kPic);
  kPic = NULL;

  if (WebPMuxGetFrame(mux, frame_number, &image, NULL,
                      &x_off, &y_off, duration) != WEBP_MUX_OK) {
    goto end;
  }
  
  output_buffer->colorspace = MODE_RGBA;
  ok = (WebPDecode(image.bytes_, image.size_, config) == VP8_STATUS_OK);

 end:
  if (!ok) {
    fprintf(stderr, "Decoding of frame #%d failed!\n", frame_number);
  } else {
    kPic = output_buffer;
  }
  return ok;
}

//------------------------------------------------------------------------------
// Main

static void Help(void) {
  printf("Usage: vwebp in_file [options]\n\n"
         "Decodes the WebP image file and visualize it using OpenGL\n"
         "Options are:\n"
         "  -version  .... print version number and exit.\n"
         "  -nofancy ..... don't use the fancy YUV420 upscaler.\n"
         "  -nofilter .... disable in-loop filtering.\n"
         "  -mt .......... use multi-threading\n"
         "  -crop <x> <y> <w> <h> ... crop output with the given rectangle\n"
         "  -scale <w> <h> .......... scale the output (*after* any cropping)\n"
         "  -h     ....... this help message.\n"
        );
}

int main(int argc, char *argv[]) {
  WebPDecoderConfig config;
  void* data;
  size_t size;
  int c;

  if (!WebPInitDecoderConfig(&config)) {
    fprintf(stderr, "Library version mismatch!\n");
    return -1;
  }

  for (c = 1; c < argc; ++c) {
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-nofancy")) {
      config.options.no_fancy_upsampling = 1;
    } else if (!strcmp(argv[c], "-nofilter")) {
      config.options.bypass_filtering = 1;
    } else if (!strcmp(argv[c], "-version")) {
      const int version = WebPGetDecoderVersion();
      printf("%d.%d.%d\n",
        (version >> 16) & 0xff, (version >> 8) & 0xff, version & 0xff);
      return 0;
    } else if (!strcmp(argv[c], "-mt")) {
      config.options.use_threads = 1;
    } else if (!strcmp(argv[c], "-crop") && c < argc - 4) {
      config.options.use_cropping = 1;
      config.options.crop_left   = strtol(argv[++c], NULL, 0);
      config.options.crop_top    = strtol(argv[++c], NULL, 0);
      config.options.crop_width  = strtol(argv[++c], NULL, 0);
      config.options.crop_height = strtol(argv[++c], NULL, 0);
    } else if (!strcmp(argv[c], "-scale") && c < argc - 2) {
      config.options.use_scaling = 1;
      config.options.scaled_width  = strtol(argv[++c], NULL, 0);
      config.options.scaled_height = strtol(argv[++c], NULL, 0);
    } else if (argv[c][0] == '-') {
      printf("Unknown option '%s'\n", argv[c]);
      Help();
      return -1;
    } else {
      file_name = argv[c];
    }
  }

  if (file_name == NULL) {
    printf("missing input file!!\n");
    Help();
    return -1;
  }
  data = ReadFile(file_name, &size);
  if (data == NULL) {
    return -1;
  }

  {
    uint32_t loop_count = 1, loop_num = 0, duration = 0;
    uint32_t flags;
    int num_images;
    WebPMuxError mux_err;
    WebPMux* const mux = WebPMuxCreate(data, size, 0, NULL);
    if (mux == NULL) {
      printf("Could not create demuxing object!\n");
      free(data);
      return -1;
    }

    mux_err = WebPMuxGetFeatures(mux, &flags);
    if (mux_err != WEBP_MUX_OK) {
      goto End;
    }

    if (flags & ~ANIMATION_FLAG) {
      fprintf(stderr, "Only extended format files containing "
                      "animation are supported for now!\n");
 End:
      free(data);
      WebPMuxDelete(mux);
      WebPFreeDecBuffer((WebPDecBuffer*)kPic);
      return -1;
    }
    mux_err = WebPMuxGetLoopCount(mux, &loop_count);
    if (mux_err != WEBP_MUX_OK && mux_err != WEBP_MUX_NOT_FOUND) {
      goto End;
    }
    num_images = 1;
    mux_err = WebPMuxNumNamedElements(mux, "image", &num_images);
    if (mux_err != WEBP_MUX_OK && mux_err != WEBP_MUX_NOT_FOUND) {
      goto End;
    }
    printf("Found %d images in file (loop_count = %d)\n",
           num_images, loop_count);

    glutInit(&argc, argv);
    for (loop_num = 0; !exiting && loop_num < loop_count; ++loop_num) {
      struct timeval now;
      uint64_t now_ms, then_ms;
      int i;

      for (i = 1; !exiting && i <= num_images; ++i) {
        if (!Decode(mux, i, &config, &duration)) {
          goto End;
        }
        gettimeofday(&now, NULL);
        now_ms = now.tv_sec * 1000 + now.tv_usec / 1000;
        then_ms = now_ms + duration;

        if (i == 1) {
          printf("Displaying [%s]: %d x %d. Press Esc to exit, 'i' for info.\n",
              file_name, kPic->width, kPic->height);
          Show(kPic);
        }
        while (now_ms < then_ms) {
          gettimeofday(&now, NULL);
          now_ms = now.tv_sec * 1000 + now.tv_usec / 1000;
          if (i > 1) {
            HandleDisplay();
          }
        }
        glutMainLoopEvent();
        sleep(0);
      }
    }
    WebPMuxDelete(mux);
  }
  if (!exiting) glutMainLoop();

  // Should only be reached when using FREEGLUT:
  free(data);
  WebPFreeDecBuffer((WebPDecBuffer*)kPic);
  return 0;
}

//------------------------------------------------------------------------------
