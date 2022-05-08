
#include "raylib.h"
#include "raymath.h"
#include <cstring>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <chrono>


struct Block
{
    Rectangle rect;
    Color color;
};


Rectangle vline2rect(Vector2 A, float height, float thickness)
{
    return Rectangle{A.x - thickness / 2,
                     A.y - thickness / 2,
                     thickness,
                     height + thickness};
}
Rectangle hline2rect(Vector2 A, float width, float thickness)
{
    return Rectangle{A.x - thickness / 2,
                     A.y - thickness / 2,
                     width + thickness,
                     thickness};
}
Rectangle vline2rect(float x, float y1, float y2, float thickness)
{
    return Rectangle{x - thickness / 2,
                     y1 - thickness / 2,
                     thickness,
                     y2 - y1 + thickness};
}
Rectangle hline2rect(float y, float x1, float x2, float thickness)
{
    return Rectangle{x1 - thickness / 2,
                     y - thickness / 2,
                     x2 - x1 + thickness,
                     thickness};
}

std::vector<Block> Box(Vector2 A, float w, float h, Color c, float thickness)
{
    return {
        { vline2rect(A.x, A.y, A.y+h, thickness), GRAY },
        { vline2rect(A.x+w, A.y, A.y+h, thickness), GRAY },
        { hline2rect(A.y, A.x, A.x+w, thickness), GRAY },
        { hline2rect(A.y+h, A.x, A.x+w, thickness), GRAY }
    };
}


////////////////////////////////////////////////////////////////////////////////


Vector2 invert(Vector2 a)
{
    return {-a.x, -a.y};
}

Vector2& operator+= (Vector2& lhs, const Vector2& rhs)
{
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    return lhs;
}

Vector2& operator-= (Vector2& lhs, const Vector2& rhs)
{
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    return lhs;
}

Vector2& operator*= (Vector2& lhs, float rhs)
{
    lhs.x *= rhs;
    lhs.y *= rhs;
    return lhs;
}

Vector2& operator/= (Vector2& lhs, float rhs)
{
    lhs.x /= rhs;
    lhs.y /= rhs;
    return lhs;
}

Vector2 operator+ (Vector2 lhs, Vector2 rhs)
{
    lhs += rhs;
    return lhs;
}

Vector2 operator- (Vector2 lhs, Vector2 rhs)
{
    lhs -= rhs;
    return lhs;
}

Vector2 operator* (Vector2 lhs, float rhs)
{
    lhs *= rhs;
    return lhs;
}

Vector2 operator/ (Vector2 lhs, float rhs)
{
    lhs /= rhs;
    return lhs;
}


float dot(Vector2 a, Vector2 b)
{
    return a.x*b.x + a.y*b.y;
}

float length2(Vector2 a)
{
    return dot(a, a);
}

float length(Vector2 a)
{
    return sqrtf(length2(a));
}

float distance2(Vector2 a, Vector2 b)
{
    return length2(a - b);
}

float distance(Vector2 a, Vector2 b)
{
    return sqrtf(distance2(a, b));
}

Vector2 norm(Vector2 a)
{
    float len = length(a);
    if (len != 0)
        return {a.x/len, a.y/len};
    else
        return {0, 0};
}

Vector2 proj(Vector2 a, Vector2 b)
{
    float s = dot(a, b)/dot(b, b);
    return {b.x*s, b.y*s};
}

std::ostream& operator<< (std::ostream& os, const Vector2& v)
{
    os << std::fixed << std::setprecision(2) << "{.x=" << v.x << ", .y=" << v.y << "}";
    return os;
}

std::string to_string(const Vector2& v)
{
    std::stringstream ss;
    ss << v;
    return ss.str();
}


////////////////////////////////////////////////////////////////////////////////


std::random_device rand_dev;
std::default_random_engine dre1(rand_dev());

float randrange(float from, float to)
{
    std::uniform_real_distribution<float> uniform_dist(from, to);
    return uniform_dist(dre1);
}

float fit_to_range(float value, std::pair<float, float> range)
{
    if (value < range.first) {
        return range.first;
    } else if (value > range.second) {
        return range.second;
    } else {
        return value;
    }
}

////////////////////////////////////////////////////////////////////////////////


struct Atom
{
    static constexpr float KM = 200;

    float r;

    Vector2 pos;
    Vector2 vel;
    Vector2 acc = {0, 0};

    void set_T(float T)
    {
        vel = norm(vel) * sqrtf(T*KM);
    }

    float get_T() const
    {
        return length2(vel)/KM;
    }

    void update_step(Vector2 global_acc, float dt)
    {
        Vector2 a = acc + global_acc;
        vel += a*dt;
        pos += vel*dt;
    }
};

std::ostream& operator<< (std::ostream& os, const Atom& a)
{
    os << "{.r=" << a.r
        << ", .pos=" << a.pos
        << ", .vel=" << a.vel
        << ", .T=" << a.get_T()
        << "}";
    return os;
}

std::string to_string(const Atom& a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}

//TODO: different impact models

void impact(Atom& atom1, Atom& atom2)
{
    float d2 = distance2(atom1.pos, atom2.pos);
    float max_distance = atom1.r + atom2.r;
    if (d2 > max_distance*max_distance) {
        return;
    }

    // vector from center of atom1 to center of atom2 goes through impact place
    Vector2 s12 = atom2.pos - atom1.pos;
    Vector2 s21 = invert(s12);

    // normal vectors to impact line - atoms will exchange velocity projections on this vector
    Vector2 n12 = norm(s12);
    Vector2 n21 = invert(n12);

    // velocity of atom1 in projection on vector from atom1 to atom2 - then goes to atom2
    Vector2 v1n12 = proj(atom1.vel, s12);
    // velocity of atom2 in projection on vector from atom2 to atom1 - then goes to atom1
    Vector2 v2n21 = proj(atom2.vel, s21);
    // left velocity - saved
    Vector2 v1p = atom1.vel - v1n12;
    Vector2 v2p = atom2.vel - v2n21;

    // exchange impulse
    atom1.vel = v1p + v2n21;
    atom2.vel = v2p + v1n12;

    float d = sqrt(d2);
    float r = max_distance/2 - d/2;
    r += 1.0;

    atom1.pos += n21 * r;
    atom2.pos += n12 * r;
}

//TODO: grid class, expandable, nailed to model coordinates

using grid_t = std::vector<std::vector<std::vector<Atom*>>>;

grid_t new_grid(size_t w, size_t h)
{
    grid_t g;
    for (size_t i = 0; i < h; i++) {
        g.push_back({});
        for (size_t j = 0; j < w; j++) {
            g.back().push_back({});
        }
    }
    return g;
}

void print_grid(const grid_t& g)
{
    printf("Grid width=%lu, height=%lu:\n", g.empty() ? 0 : g[0].size(), g.size());
    size_t iy = 0;
    for (const auto& row : g) {
        size_t ix = 0;
        for (const auto& cell : row) {
            if (cell.size() > 0) {
                printf("    Cell %lu %lu:\n", iy, ix);
                for (const Atom* atom : cell) {
                    printf("        %p: %s\n", atom, to_string(*atom).c_str());
                }
            }
            ix++;
        }
        iy++;
    }
}

struct MyModel
{
    float time_scale = 10.0;

    static constexpr float T_MAX = 1000;
    static constexpr float T_MIN = 273.15;

    float R_atom = 1000;

    Vector2 g_acc = {0, -9.81 / 5};

    std::vector<Atom> atoms;

    float width = 400 * R_atom;
    float height = 400 * R_atom;

    grid_t grid;

    size_t grid_row_size = 100;
    size_t grid_col_size = 100;

    float grid_cell_width() const
    {
        return width / grid_row_size;
    }

    float grid_cell_height() const
    {
        return height / grid_col_size;
    }

    float max_atoms_in_cell() const
    {
        return grid_cell_height() * grid_cell_height() / (4*R_atom*R_atom);
    }

    Atom& create_atom(Vector2 pos, Vector2 v)
    {
        Atom atom{R_atom, pos, v};

        atoms.push_back(atom);
        Atom& a = atoms.back();

        index_atom(grid, &a);

        return a;
    }

    void index_atom(grid_t& g, Atom* atom)
    {
        size_t ix = atom->pos.x / grid_cell_width();
        size_t iy = atom->pos.y / grid_cell_height();
        g[iy][ix].push_back(atom);
    }

    void print_atoms()
    {
        printf("atoms (count=%lu):\n", atoms.size());
        for (const auto& atom : atoms) {
            printf("    %p: %s\n", &atom, to_string(atom).c_str());
        }
    }

    Atom& create_test_atom(float ix, float iy, float dx, float dy)
    {
        float x = (ix + 0.5)*R_atom*2;
        float y = (iy + 0.5)*R_atom*2;
        x = fit_to_range(x, {R_atom, width - R_atom});
        y = fit_to_range(y, {R_atom, height - R_atom});
        return create_atom({x, y}, {dx, dy});
    }

    void setup_debug()
    {
        time_scale = 20;

        width = R_atom * 16;
        height = R_atom * 16;

        grid_row_size = 5;
        grid_col_size = 5;

        grid = new_grid(grid_row_size, grid_col_size);

        g_acc = {0, 0};

        atoms.reserve(16);

        create_test_atom(0, 2, 500, 0);
        create_test_atom(2, 2, 0, 0);
        create_test_atom(3, 2, 0, 0);
        create_test_atom(2, 6, 0, -20);

        print_atoms();
        print_grid(grid);
    }

    void setup_random(size_t count)
    {
        grid = new_grid(grid_row_size, grid_col_size);

        atoms.reserve(count);

        for (size_t i = 0; i < count; i++) {
            float x = randrange(R_atom, width - R_atom);
            float y = randrange(R_atom, height - R_atom);
            float qx = randrange(-10, 10);
            float qy = randrange(-10, 10);
            float t = randrange(0.0, T_MIN);
            create_atom({x, y}, {qx, qy}).set_T(t);
        }
    }

    void update(float dt_seconds)
    {
        float dt = dt_seconds * time_scale;

        //TODO: max speed not exceed 1/2*R_atom

        grid_t result_grid = new_grid(grid_row_size, grid_col_size);

        int iy = 0;
        for (const auto& row : grid) {
            int ix = 0;
            for (const auto& cell : row) {
                for (Atom* atom : cell) {

                    for (int pix = ix-1; pix < ix+1; pix++) {
                        for (int piy = iy-1; piy < iy+1; piy++) {
                            if (pix >= 0 && piy >= 0 && pix < int(row.size()) && piy < int(grid.size())) {
                                for (Atom* a : grid[piy][pix]) {
                                    if (a != atom) {
                                        impact(*atom, *a);
                                    }
                                }
                            }
                        }
                    }

                    atom->update_step(g_acc, dt);

                    // borders
                    if (atom->pos.x < atom->r && atom->vel.x < 0) {
                        atom->pos.x = atom->r;
                        atom->vel.x = -atom->vel.x;
                    }
                    if (atom->pos.x > width - atom->r && atom->vel.x > 0) {
                        atom->pos.x = width - atom->r;
                        atom->vel.x = -atom->vel.x;
                    }
                    if (atom->pos.y < atom->r && atom->vel.y < 0) {
                        atom->pos.y = atom->r;
                        atom->vel.y = -atom->vel.y;
                        atom->set_T(T_MAX);
                    }
                    if (atom->pos.y > height - atom->r && atom->vel.y > 0) {
                        atom->pos.y = height - atom->r;
                        atom->vel.y = -atom->vel.y;
                        atom->set_T(T_MIN);
                    }

                    index_atom(result_grid, atom);
                }
                ix++;
            }
            iy++;
        }

        grid = result_grid;
    }

    struct GridCellDensity
    {
        Rectangle r;
        float density;
    };

    std::vector<GridCellDensity> get_cells() const
    {
        std::vector<GridCellDensity> d;
        d.reserve(grid_col_size * grid_row_size);

        size_t max_count = max_atoms_in_cell() + 1;
        float cell_h = grid_cell_height();
        float cell_w = grid_cell_width();

        int iy = 0;
        for (const auto& row : grid) {
            int ix = 0;
            for (const auto& cell : row) {

                d.push_back({
                    {
                        ix * cell_w,
                        iy * cell_h,
                        cell_w,
                        cell_h
                    },
                    float(cell.size()) / max_count
                });

                ix++;
            }
            iy++;
        }

        return d;
    }

}; // MyModel



struct MyRenderer
{
    float scale_factor = 0;

    float gradient_radius = 0;
    float gradient_scale = 1.0;

    float window_width = 0;
    float window_height = 0;

    const uint8_t gradient_center_opacity = 100;
    const uint8_t gradient_edge_opacity = 0;

    bool show_density = false;
    bool show_atoms = true;

    void setup(const MyModel& m, float screen_width, float screen_height)
    {
        scale_factor = std::max(m.height / screen_height, m.width / screen_width);

        gradient_radius = m.R_atom / scale_factor * gradient_scale;
        window_width = m.width / scale_factor;
        window_height = m.height / scale_factor;
    }

    Vector2 atom_pos_to_screen(Vector2 p)
    {
        float x = p.x/scale_factor - window_width/2 - gradient_radius;
        float y = -p.y/scale_factor + window_height/2 - gradient_radius;
        return {x, y};
    }

    uint8_t temperature_to_sprite_index(float t)
    {
        float t_min = 0;
        float t_max = MyModel::T_MAX * 2;

        if (t < t_min) {
            t = t_min;
        }
        if (t > t_max) {
            t = t_max;
        }

        int c = (t-t_min)/(t_max-t_min)*255;

        if (c < 0) {
            return 0;
        }
        if (c > 255) {
            return 255;
        }

        return uint8_t(c);
    }

    Color gradient_inner_color(uint8_t i)
    {
        return {uint8_t(i), 0, uint8_t(255-i), gradient_center_opacity};
    }

    Color gradient_outer_color(uint8_t i)
    {
        return {uint8_t(i), 0, uint8_t(255-i), gradient_edge_opacity};
    }

    std::vector<Texture2D> atom_colors;

    void prepare()
    {
        for (size_t i = 0; i < 256; i++) {
            Image radial_gradient = GenImageGradientRadial(gradient_radius * 2,
                                                        gradient_radius * 2,
                                                        0.0f,
                                                        gradient_inner_color(i),
                                                        gradient_outer_color(i));
            atom_colors.push_back(LoadTextureFromImage(radial_gradient));
            UnloadImage(radial_gradient);
        }
    }

    Rectangle from_model(Rectangle r)
    {
        return {
            r.x / scale_factor - window_width/2,
            -r.y / scale_factor + window_height/2 - r.height/scale_factor,
            r.width / scale_factor,
            r.height / scale_factor
        };
    }
    Color cell_color(float density)
    {
        int c = 40 + int(80 * density);
        uint8_t u = c <= 255 ? 255-c : 0;
        return {u, u, u, 255};
    }

    void render(const MyModel& m)
    {
        if (show_density) {
            auto cells = m.get_cells();

            for (const auto& cell : cells) {
                Color c = cell_color(cell.density);
                auto r = from_model(cell.r);

                DrawRectangleRec(r, c);
            }
        } else {
            DrawRectangleRec({-window_width/2, -window_height/2, window_width, window_height}, LIGHTGRAY);
        }

        if (show_atoms) {
            for (const auto& atom : m.atoms) {
                float t = atom.get_T();
                uint8_t color_index = temperature_to_sprite_index(t);
                Texture2D& tx = atom_colors[color_index];
                Vector2 pos = atom_pos_to_screen(atom.pos);
                DrawTextureEx(tx, pos, 0, 1.0, WHITE);
            }
        }
    }

}; // MyRenderer



using Clock = std::chrono::steady_clock;


int main(void)
{
    MyModel model;

    //model.setup_debug();
    model.setup_random(10000);

    MyRenderer r;
    r.show_density = false;
    r.show_atoms = true;
    r.gradient_scale = 4.0;


    int screen_width = 800;
    int screen_height = 600;
    int max_fps = 120;

    InitWindow(screen_width, screen_height, "2D atomic gas");

    r.setup(model, screen_width, screen_height);
    r.prepare();

    Camera2D camera = { {0} };
    camera.target = {0, 0};
    camera.offset = (Vector2){ screen_width / 2.0f, screen_height / 2.0f };
    camera.rotation = 0.0f;
    camera.zoom = 1.0f;

    SetTargetFPS(max_fps);

    Clock::time_point last_ts = Clock::now();

    while (!WindowShouldClose()) {
        Clock::time_point now = Clock::now();
        std::chrono::duration<float> dt = now - last_ts;
        last_ts = now;

        model.update(dt.count());

        BeginDrawing();
        ClearBackground(BLACK);
        BeginMode2D(camera);

        r.render(model);

        EndMode2D();
        DrawFPS(20, 5);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
