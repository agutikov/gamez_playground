
#include <cstdio>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <chrono>
#include <limits>
#include <exception>
#include <iostream>

#include <fmt/core.h>

#include "raylib.h"
#include "raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#include "lib/geom.h"
#include "lib/vector.h"
#include "lib/random.h"


using namespace std::chrono_literals;

using Clock = std::chrono::steady_clock;


struct Atom
{
    static constexpr float K = 30000;

    float r;

    Vector2 pos;
    Vector2 vel;
    Vector2 acc = {0, 0};

    const float m = 1.0;

    void set_T(float T)
    {
        vel = norm(vel) * sqrtf(T*K*m);
    }

    float get_T() const
    {
        return length2(vel)/(K*m);
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
//TODO: mass
//TODO: different R

void impact(Atom& atom1, Atom& atom2)
{
    float d2 = distance2(atom1.pos, atom2.pos);
    float max_distance = atom1.r + atom2.r;
    if (d2 > max_distance*max_distance) {
        return;
    }

    // vector from center of atom1 to center of atom2 goes through impact place
    Vector2 s12 = atom2.pos - atom1.pos;
    Vector2 s21 = -s12;

    // normal vectors to impact line - atoms will exchange velocity projections on this vector
    Vector2 n12 = norm(s12);
    Vector2 n21 = -n12;

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

    bool top_exchange_temperature = true;
    bool bottom_exchange_temperature = true;
    float T_top = 10;
    float T_bottom = 2000;

    float R_atom = 1000;

    Vector2 g_acc = {0, -9.81};

    std::vector<Atom> atoms;

    float width;
    float height;

    grid_t grid;
    std::vector<std::vector<Rectangle>> grid_cells;

    float grid_cell_R_count = 4.0;
    size_t grid_row_size;
    size_t grid_col_size;

    float max_atoms_in_cell;

    stat_value<float>speed_max_values;
    stat_value<float>speed_mean_values;
    stat_value<float>speed_min_values;

    decltype(speed_max_values)::stats_t speed_max = {};
    decltype(speed_mean_values)::stats_t speed_mean = {};
    decltype(speed_min_values)::stats_t speed_min = {};

    float last_T_max = 0;
    float last_speed_max = 0;
    float last_R_min = 0;

    Clock::time_point stats_timestamp = Clock::now();
    std::chrono::milliseconds stats_period = 2s;

    size_t max_atoms_count = 10000;
    float new_atoms_rate = 0;
    float accumulated_time = 0;

    float grid_cell_width() const
    {
        return width / grid_row_size;
    }

    float grid_cell_height() const
    {
        return height / grid_col_size;
    }

    float grid_cell_density(size_t x, size_t y) const
    {
        return grid[y][x].size() / max_atoms_in_cell;
    }

    Atom& create_atom(Vector2 pos, Vector2 v)
    {
        if (atoms.size() >= max_atoms_count) {
            throw std::runtime_error("Exceeds max number of atoms");
        }

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

    MyModel(Vector2 dimensions)
    {
        width = dimensions.x * R_atom;
        height = dimensions.y * R_atom;

        grid_row_size = width / (R_atom * grid_cell_R_count);
        grid_col_size = height / (R_atom * grid_cell_R_count);
    }

    void basic_setup()
    {
        atoms.reserve(max_atoms_count);

        grid = new_grid(grid_row_size, grid_col_size);
        setup_grid_cells();
        max_atoms_in_cell = grid_cell_height() * grid_cell_height() / (4*R_atom*R_atom);
    }

    void setup_debug()
    {
        time_scale = 20;

        width = R_atom * 16;
        height = R_atom * 16;

        grid_row_size = 5;
        grid_col_size = 5;

        g_acc = {0, 0};

        basic_setup();

        create_test_atom(0, 2, 500, 0);
        create_test_atom(2, 2, 0, 0);
        create_test_atom(3, 2, 0, 0);
        create_test_atom(2, 6, 0, -20);

        print_atoms();
        print_grid(grid);
    }

    void generate_random(size_t count)
    {
        for (size_t i = 0; i < count; i++) {
            float x = randrange(R_atom, width - R_atom);
            float y = randrange(R_atom, height - R_atom);
            float qx = randrange(-10, 10);
            float qy = randrange(-10, 10);
            float t = randrange(0.0, T_top);
            create_atom({x, y}, {qx, qy}).set_T(t);
        }
    }

    void update_stats()
    {
        auto t = Clock::now();
        if (t < stats_timestamp + stats_period) {
            return;
        }
        stats_timestamp = t;

        speed_max = speed_max_values.produce();
        speed_min = speed_min_values.produce();
        speed_mean = speed_mean_values.produce();
    }

    void update_new_atoms(float dt)
    {
        if (atoms.size() >= max_atoms_count) {
            return;
        }
        if (new_atoms_rate <= 0) {
            return;
        }
        accumulated_time += dt;
        size_t count = new_atoms_rate * accumulated_time;
        if (count <= 0) {
            return;
        }
        accumulated_time -= count / new_atoms_rate;

        while (count--) {
            create_atom({width - 2*R_atom, height - 10*R_atom}, {-4000, 0});
        }
    }

    void update_time_scale(float dt_seconds)
    {
        float max_move = last_R_min / 2;
        float dt = dt_seconds * time_scale;
        float next_move = last_speed_max * dt;
        if (next_move > max_move) {
            // max_move == next_move = last_speed_max * dt_seconds * time_scale
            time_scale = max_move / (last_speed_max * dt_seconds);
        }
    }

    void update(float dt_seconds)
    {
        stat_value<float>temperature;
        stat_value<float>speed;
        stat_value<float>radius; // in case if radius will be variable

        // max speed must not exceed (min(R)/2)/dt, adjust time_scale instead
        update_time_scale(dt_seconds);
        float dt = dt_seconds * time_scale;

        update_new_atoms(dt);

        grid_t result_grid = new_grid(grid_row_size, grid_col_size);

        int iy = 0;
        for (const auto& row : grid) {
            int ix = 0;
            for (const auto& cell : row) {
                for (Atom* atom : cell) {

                    //TODO: rewrite and optimize

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
                        if (bottom_exchange_temperature) {
                            atom->set_T(T_bottom);
                        }
                    }
                    if (atom->pos.y > height - atom->r && atom->vel.y > 0) {
                        atom->pos.y = height - atom->r;
                        atom->vel.y = -atom->vel.y;
                        if (top_exchange_temperature) {
                            atom->set_T(T_top);
                        }
                    }

                    temperature.consume(atom->get_T());
                    speed.consume(length(atom->vel));
                    radius.consume(atom->r);

                    index_atom(result_grid, atom);
                }
                ix++;
            }
            iy++;
        }

        grid = result_grid;

        auto speed_stats = speed.produce();

        last_T_max = temperature.produce().max;
        last_speed_max = speed_stats.max;
        last_R_min = radius.produce().min;

        speed_min_values.consume(speed_stats.min / R_atom);
        speed_mean_values.consume(speed_stats.mean / R_atom);
        speed_max_values.consume(speed_stats.max / R_atom);

        update_stats();
    }

    void setup_grid_cells()
    {
        float cell_h = grid_cell_height();
        float cell_w = grid_cell_width();

        for (size_t y = 0; y < grid_col_size; y++) {
            grid_cells.push_back({});
            for (size_t x = 0; x < grid_row_size; x++) {
                grid_cells.back().push_back({
                    x * cell_w,
                    y * cell_h,
                    cell_w,
                    cell_h
                });
            }
        }
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

    bool draw_density_map = false;
    bool draw_atoms = true;

    std::vector<Texture2D> atom_colors;

    Texture2D selected_atom;
    Color selected_color = BLACK;

    float margin = 10;

    std::vector<std::vector<Rectangle>> grid_cells;

    Rectangle background;

    void setup(const MyModel& m, float screen_width, float screen_height)
    {
        scale_factor = std::min(screen_height / m.height, screen_width / m.width);

        gradient_radius = m.R_atom * scale_factor * gradient_scale;
        window_width = m.width * scale_factor;
        window_height = m.height * scale_factor;

        setup_grid_cells(m);

        background = {0, 0, window_width, window_height};
    }

    Vector2 atom_pos_to_screen(Vector2 p)
    {
        float x = p.x * scale_factor - gradient_radius;
        float y = -p.y * scale_factor + window_height - gradient_radius;
        return {x, y};
    }

    Color gradient_inner_color(uint8_t i)
    {
        return {uint8_t(i), 0, uint8_t(255-i), gradient_center_opacity};
    }

    Color gradient_outer_color(uint8_t i)
    {
        return {uint8_t(i), 0, uint8_t(255-i), gradient_edge_opacity};
    }

    Texture2D create_radial_gradient(Color inner, Color outer, float fat)
    {
        Image radial_gradient = GenImageGradientRadial(gradient_radius * 2,
                                                    gradient_radius * 2,
                                                    fat,
                                                    inner,
                                                    outer);
        Texture2D txt = LoadTextureFromImage(radial_gradient);
        UnloadImage(radial_gradient);
        return txt;
    }

    void prepare()
    {
        for (size_t i = 0; i < 256; i++) {
            atom_colors.push_back(
                create_radial_gradient(
                    gradient_inner_color(i),
                    gradient_outer_color(i),
                    0.0f
                )
            );
        }
        selected_atom = create_radial_gradient(selected_color, BLANK, 0.0f);
    }

    void free()
    {
        for (auto& txt : atom_colors) {
            UnloadTexture(txt);
        }
        UnloadTexture(selected_atom);
    }

    Rectangle to_screen(Rectangle r)
    {
        return {
            r.x * scale_factor,
            window_height - (r.y + r.height) * scale_factor,
            r.width * scale_factor,
            r.height * scale_factor
        };
    }

    void setup_grid_cells(const MyModel& m)
    {
        for (size_t y = 0; y < m.grid_col_size; y++) {
            grid_cells.push_back({});
            for (size_t x = 0; x < m.grid_row_size; x++) {
                grid_cells.back().push_back(to_screen(m.grid_cells[y][x]));
            }
        }
    }

    Color cell_color(float density)
    {
        int c = 40 + int(80 * density);
        uint8_t u = c <= 255 ? 255-c : 0;
        return {u, u, u, 255};
    }

    void render(const MyModel& m)
    {
        if (draw_density_map) {
            for (size_t y = 0; y < m.grid_col_size; y++) {
                for (size_t x = 0; x < m.grid_row_size; x++) {
                    Color c = cell_color(m.grid_cell_density(x, y));

                    DrawRectangleRec(grid_cells[y][x], c);
                }
            }
        } else {
            DrawRectangleRec(background, WHITE);
        }

        if (draw_atoms) {
            for (auto it = m.atoms.begin() + 1; it < m.atoms.end(); ++it) {
                const auto& atom = *it;

                float t = atom.get_T() / m.last_T_max;
                t = fit_to_range(t, {0.0, 1.0});
                uint8_t color_index = t * 255;

                const auto& tx = atom_colors[color_index];
                Vector2 pos = atom_pos_to_screen(atom.pos);

                DrawTextureEx(tx, pos, 0, 1.0, WHITE);
            }

            Vector2 pos = atom_pos_to_screen(m.atoms[0].pos);
            DrawTextureEx(selected_atom, pos, 0, 1.0, WHITE);
        }
    }

    void render_text(const MyModel& m)
    {
        DrawText(fmt::format("N={}", m.atoms.size()).c_str(), 5, 40, 12, GRAY);
        DrawText(fmt::format("Time scale={:.3f}", m.time_scale).c_str(), 5, 60, 12, GRAY);
        DrawText("Speed/R_atom:", 5, 80, 12, GRAY);
        DrawText(("  min=" + std::string(m.speed_min)).c_str(), 5, 100, 12, GRAY);
        DrawText(("  mean=" + std::string(m.speed_mean)).c_str(), 5, 120, 12, GRAY);
        DrawText(("  max=" + std::string(m.speed_max)).c_str(), 5, 140, 12, GRAY);
    }

}; // MyRenderer


//TODO: mark one or several atom with color
//TODO: record coordinates track of some atom
//TODO: random generator

//TODO: imgui

//TODO: cylinder and tor surface

//TODO: window resize and fullscreen
//TODO: zoom with zoomed gradients


int main(void)
{
    MyModel model({100, 100});

    model.basic_setup();

    model.top_exchange_temperature = false;
    model.bottom_exchange_temperature = false;
    model.max_atoms_count = 1000;

    //model.setup_debug();
    model.generate_random(1000);
    //model.g_acc = {0, -9.81 * 40};
    //model.new_atoms_rate = 1.0;

    MyRenderer r;
    r.draw_density_map = false;
    r.draw_atoms = true;
    r.gradient_scale = 1.0;

    int screen_width = 900;
    int screen_height = 700;
    int max_fps = 120;
    IntRectangle viewport{200, 10, screen_width-210, screen_height-20};
    Rectangle view = viewport;
    Vector2 view_offset{ view.x, view.y };

    InitWindow(screen_width, screen_height, "2D Particles");
    SetExitKey(0);

    r.setup(model, viewport.width, viewport.height);
    r.prepare();

    Camera2D camera = { {0} };
    camera.target = {0, 0};
    camera.offset = view_offset;
    camera.rotation = 0.0f;
    camera.zoom = 1.0f;

    SetTargetFPS(max_fps);

    Clock::time_point last_ts = Clock::now() - 10ms;

    Vector2 mouse_drag_from = {0, 0};
    bool drag_on = false;

    bool play = true;
    Rectangle playbutton_rec{16, 256, 32, 32};

    while (!WindowShouldClose()) {
        Clock::time_point now = Clock::now();
        std::chrono::duration<float> dt = now - last_ts;
        last_ts = now;

        if (play) {
            model.update(dt.count());
        }

        float wheel = GetMouseWheelMove();
        Vector2 mouse = GetMousePosition();
        if (wheel != 0.0 && CheckCollisionPointRec(mouse, view)) {
            Vector2 target_mouse = GetScreenToWorld2D(mouse, camera);
            camera.target = target_mouse;
            camera.offset = mouse;
            camera.zoom += wheel * 0.05;
        }
        if (drag_on) {
            Vector2 move = mouse - mouse_drag_from;
            camera.offset += move;
            mouse_drag_from = mouse;
        }
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON) && !drag_on) {
            mouse_drag_from = mouse;
            drag_on = true;
        }
        if (!IsMouseButtonDown(MOUSE_LEFT_BUTTON) && drag_on) {
            drag_on = false;
        }

        if (IsKeyPressed(KEY_SPACE)) {
            play = !play;
        }

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode2D(camera);
        BeginScissorMode(viewport.x, viewport.y, viewport.width, viewport.height);
        ClearBackground(DARKGRAY);

        r.render(model);

        EndMode2D();
        EndScissorMode();

        guiIconName play_icon = play
            ? guiIconName::RAYGUI_ICON_PLAYER_PAUSE
            : guiIconName::RAYGUI_ICON_PLAYER_PLAY;
        if (GuiButton(playbutton_rec, GuiIconText(play_icon, ""))) {
            printf("%s\n", play ? "Pause" : "Play");
            play = !play;
        }


        DrawFPS(20, 5);
        r.render_text(model);

        EndDrawing();
    }

    r.free();

    CloseWindow();

    return 0;
}
