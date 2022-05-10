#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>
#include <exception>
#include <iostream>

#include "raylib.h"
#include "raymath.h"

#include "PerlinNoise.hpp"

#include "lib/vector.h"
#include "lib/random.h"
#include "lib/geom.h"



struct far_stars_cfg
{
    bool enable;
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
    uint64_t seed;
    perlin_cfg perlin;
    size_t min_count;
    size_t max_count;
};

struct close_stars_cfg
{
    bool enable;
};

struct space_cfg
{
    far_stars_cfg far_stars;
    nebulae_cfg nebulae;
    close_stars_cfg close_stars;
};


void generate_far_stars(Image* img, far_stars_cfg cfg)
{
    if (!cfg.enable) {
        return;
    }

}


int randrange(int from, int to)
{
    static std::random_device rand_dev;
    static std::default_random_engine dre(rand_dev());

    std::uniform_int_distribution<int> uniform_dist(from, to);
    return uniform_dist(dre);
}

Color rand_color()
{
    uint32_t value = randrange(0, (0x1 << 24) - 1);
    return {
        uint8_t(value),
        uint8_t(value >> 8),
        uint8_t(value >> 16),
        255
    };
}

Vector2 rand_vec(float x_max, float y_max)
{
    return {
        randrange(0.0f, x_max),
        randrange(0.0f, y_max)
    };
}

Rectangle ImageRec(const Image& img)
{
    return { 0, 0, float(img.width), float(img.height)};
}


void image_generate_perlin(Image* img, perlin_cfg cfg)
{
    const siv::PerlinNoise perlin{cfg.seed};

    const double fx = (cfg.frequency / img->width);
    const double fy = (cfg.frequency / img->height);

    for (int y = 0; y < img->height; y++) {
        for (int x = 0; x < img->width; x++) {
            float noise = perlin.octave2D_01((x * fx), (y * fy), cfg.octaves, cfg.persistence);
            Color c = Fade(WHITE, noise);
            ImageDrawPixel(img, x, y, c);
        }
    }
}

//TODO: RenderTexture2D and BeginBlendMode

//TODO: cache perlin


void generate_nebulae(Image* img, nebulae_cfg cfg)
{
    if (!cfg.enable) {
        return;
    }

    Image perlin = GenImageColor(img->width*2, img->height*2, BLANK);
    image_generate_perlin(&perlin, cfg.perlin);

    size_t count = randrange(int(cfg.min_count), int(cfg.max_count));

    while (count--) {
        float density = randrange(0.01f, 0.2f);
        Color c = rand_color();
        c = Fade(c, density);

        float scale = 1 << randrange(int(1), int(5));

        fmt::print("\ndensity={}, scale={}\n", density, scale);

        float sample_width = img->width / scale;
        float sample_height = img->height / scale;
        Vector2 offset = rand_vec(
            float(perlin.width - sample_width),
            float(perlin.height - sample_height)
        );
        Rectangle perlin_sample_rec = {offset.x, offset.y, sample_width, sample_height};

        Image perlin_sample = ImageFromImage(perlin, perlin_sample_rec);

        std::cout << "img=" << ImageRec(*img) << std::endl;
        std::cout << "perlin=" << ImageRec(perlin) << std::endl;
        std::cout << "perlin_sample_rec=" << perlin_sample_rec << std::endl;

        ImageResize(&perlin_sample, img->width, img->height);

        std::cout << "perlin_sample=" << ImageRec(perlin_sample) << std::endl;

        ImageDraw(img, perlin_sample, ImageRec(perlin_sample), ImageRec(*img), c);

        UnloadImage(perlin_sample);
    }

    UnloadImage(perlin);
}

void generate_close_stars(Image* img, close_stars_cfg cfg)
{
    if (!cfg.enable) {
        return;
    }
}

Texture2D generate_space(int width, int height, space_cfg cfg)
{
    Image img = GenImageColor(width, height, BLANK);

    generate_far_stars(&img, cfg.far_stars);
    generate_nebulae(&img, cfg.nebulae);
    generate_close_stars(&img, cfg.close_stars);

    ImageDrawRectangleLines(&img, {100, 100, 100, 100}, 10, GREEN);

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


int main(void)
{
    int screen_width = 600;
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

    int background_width = screen_width;
    int background_height = screen_height;

    space_cfg cfg = {
        .far_stars = {
            .enable = true
        },
        .nebulae = {
            .enable = true,
            .seed = 0,
            .perlin = {
                .seed = 0,
                .octaves = 8,
                .frequency = 32.0,
                .persistence = 0.5
            },
            .min_count = 1,
            .max_count = 6
        },
        .close_stars = {
            .enable = true
        }
    };

    Texture2D background  = generate_space(background_width, background_height, cfg);

    while (!WindowShouldClose()) {
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
