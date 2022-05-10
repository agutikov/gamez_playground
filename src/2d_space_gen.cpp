#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>
#include <exception>
#include <iostream>
#include <cmath>

#include "raylib.h"
#include "raymath.h"

#include "PerlinNoise.hpp"

#include "lib/vector.h"
#include "lib/random.h"
#include "lib/geom.h"



struct far_stars_cfg
{
    bool enable;
    float brightness;
    float density;
};

struct perlin_cfg
{
    siv::PerlinNoise::seed_type seed;
    std::int32_t octaves;
    double frequency;
    float persistence;
};

struct nebulae_cfg
{
    bool enable;
    perlin_cfg perlin;
    size_t count;
    float density;
    float falloff;
};

struct close_stars_cfg
{
    bool enable;
    float density;
    float core_falloff;
    float halo_falloff;
};

struct space_cfg
{
    uint64_t seed;
    far_stars_cfg far_stars;
    nebulae_cfg nebulae;
    close_stars_cfg close_stars;
};


void generate_far_stars(Image* img, far_stars_cfg cfg)
{
    if (!cfg.enable) {
        return;
    }

    size_t count = img->width * img->height * cfg.density;

    while (count--) {
        Vector2 pos = rand_vec(img->width, img->height);

        float v = std::log(1 - randrange(0.0000001, 0.9999999)) * -cfg.brightness;
        Color c = ColorFromNormalized(Vector4{v, v, v, 1.0});

        ImageDrawPixelV(img, pos, c);
    }
}


struct perlin_2d
{
    perlin_2d(perlin_cfg cfg, float scale) :
        cfg(cfg),
        perlin(cfg.seed),
        scale(scale)
    {}

    float operator()(Vector2 p) const
    {
        return perlin.octave2D_01(p.x * scale, p.y * scale, cfg.octaves, cfg.persistence);
    }

    const perlin_cfg cfg;
    const siv::PerlinNoise perlin;
    const float scale;
};

void image_generate_perlin(Image* img, perlin_cfg cfg)
{
    //const float scale = std::min((cfg.frequency / img->width), (cfg.frequency / img->height));
    const float scale = 1.0f / 200;
    perlin_2d perlin(cfg, scale);

    for (int y = 0; y < img->height; y++) {
        for (int x = 0; x < img->width; x++) {
            Vector2 p{float(x), float(y)};
            float noise = perlin(p);
            //noise = pow(noise + 0.2, 5);
            Color c = Fade(WHITE, noise);
            ImageDrawPixel(img, x, y, c);
        }
    }
}



void generate_nebulae(Image* img, nebulae_cfg cfg)
{
    if (!cfg.enable) {
        return;
    }

    size_t count = randrangeint(1, cfg.count);

    while (count--) {
        Image tmp = GenImageColor(img->width, img->height, BLANK);

        Color color = rand_color();
        Vector2 offset = rand_vec(img->width*100, img->height*100);

        float scale = 1.0f / pow(2.0f, randrange(7.0f, 11.0f));

        perlin_2d perlin(cfg.perlin, scale);

        for (int y = 0; y < tmp.height; y++) {
            for (int x = 0; x < tmp.width; x++) {
                Vector2 p{float(x), float(y)};
                float noise = perlin(p + offset);
                noise = pow(noise + cfg.density, cfg.falloff);
                Color c = Fade(color, noise);
                ImageDrawPixel(&tmp, x, y, c);
            }
        }

        ImageDraw(img, tmp, ImageRec(tmp), ImageRec(*img), WHITE);

        UnloadImage(tmp);
    }
}



void generate_close_stars(Image* img, close_stars_cfg cfg)
{
    if (!cfg.enable) {
        return;
    }

    Image tmp = GenImageColor(img->width, img->height, BLANK);

    std::vector<float> R;
    std::vector<std::pair<Color, Color>> colors;
    std::vector<Vector2> centers;

    size_t count = img->width * img->height * cfg.density;
    count = randrangeint(count/2, count);

    for (size_t i = 0; i < count; i++) {
        Color core_color = ColorAlphaBlend(WHITE, Fade(rand_color(), 0.3), WHITE);
        Color halo_color = ColorAlphaBlend(WHITE, Fade(rand_color(), 0.4), WHITE);
        colors.push_back({core_color, halo_color});
        centers.push_back(rand_vec(img->width, img->height));
        R.push_back(randrange(0.1, 4.0));
    }

    for (int y = 0; y < tmp.height; y++) {
        for (int x = 0; x < tmp.width; x++) {
            Color c = BLANK;
            Vector2 p{float(x), float(y)};

            for (size_t i = 0; i < count; i++) {
                float R_core = R[i];

                auto [core_color, halo_color] = colors[i];
                Vector2 center = centers[i];

                float d = length(p - center);

                if (d > R_core) {
                    float e1 = exp(-(d - R_core)*cfg.core_falloff);
                    float e2 = exp(-(d - R_core)*cfg.halo_falloff);

                    Color c1 = Fade(core_color, e1);
                    Color c2 = Fade(halo_color, e2);
                    Color c3 = ColorAlphaBlend(c1, c2, WHITE);
                    c = ColorAlphaBlend(c, c3, WHITE);
                } else {
                    c = core_color;
                }
            }

            ImageDrawPixel(&tmp, x, y, c);
        }
    }

    ImageDraw(img, tmp, ImageRec(tmp), ImageRec(*img), WHITE);

    UnloadImage(tmp);
}

Texture2D generate_space(int width, int height, space_cfg cfg)
{
    Image img = GenImageColor(width, height, BLANK);

    generate_far_stars(&img, cfg.far_stars);
    generate_nebulae(&img, cfg.nebulae);
    generate_close_stars(&img, cfg.close_stars);

    //image_generate_perlin(&img, cfg.nebulae.perlin);

    //ImageDrawRectangleLines(&img, {100, 100, 100, 100}, 10, GREEN);

    Texture2D txt = LoadTextureFromImage(img);
    UnloadImage(img);
    return txt;
}

//TODO: compare different noise generators
// https://github.com/Reputeless/PerlinNoise
// https://github.com/hrhwilliams/noisy
// https://github.com/sol-prog/Perlin_Noise

//TODO: GPU perlin and texture generation
// https://wwwtyro.net/2016/10/22/2D-space-scene-procgen.html
// https://github.com/wwwtyro/space-2d/blob/gh-pages/nebula.js


//TODO: scale, drag, rotate
//TODO: window resize, maximize and back and fullscreen
//TODO: UI

//TODO: generate local nebulae
//      - mask with perlin + circle gradient, or ellipse gradient (resized)


//TODO: make generated image 100% defined by config


int main(void)
{
    int screen_width = 800;
    int screen_height = 600;
    int max_fps = 30;

    SetTraceLogLevel(LOG_WARNING);

    InitWindow(screen_width, screen_height, "2D generator");

    SetTargetFPS(max_fps);

    Camera2D camera = { {0} };
    camera.target = {0, 0};
    camera.offset = {0, 0};
    camera.rotation = 0.0f;
    camera.zoom = 1.0f;

    int background_width = 1024;
    int background_height = 1024;

    space_cfg cfg = {
        .seed = 0,
        .far_stars = {
            .enable = true,
            .brightness = 0.3,
            .density = 0.001
        },
        .nebulae = {
            .enable = true,
            .perlin = {
                .seed = 0,
                .octaves = 8,
                .frequency = 2.0,
                .persistence = 0.5
            },
            .count = 4,
            .density = 0.1,
            .falloff = 6,
        },
        .close_stars = {
            .enable = true,
            .density = 1.0f/100000,
            .core_falloff = 0.9,
            .halo_falloff = 0.35,
        }
    };

    Texture2D background  = generate_space(background_width, background_height, cfg);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_F11)) {
            int display = GetCurrentMonitor();
            if (IsWindowFullscreen()) {
                SetWindowSize(screen_width, screen_height);
                EnableCursor();
            } else {
                SetWindowSize(GetMonitorWidth(display), GetMonitorHeight(display));
                DisableCursor();
            }
            ToggleFullscreen();
        }

        if (IsKeyPressed(KEY_SPACE) | IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {
            UnloadTexture(background);
            cfg.nebulae.perlin.seed++;
            background = generate_space(background_width, background_height, cfg);
        }

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode2D(camera);

        DrawTextureEx(background, {0, 0}, 0.0, 1.0, WHITE);

        EndMode2D();

        DrawFPS(20, 5);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}
