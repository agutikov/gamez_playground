
#include "raylib.h"
#include "raymath.h"
#include <cstring>
#include <cstdio>
#include <vector>
#include <cmath>
#include <random>
#include <string>


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


float gradient_radius = 10;
uint8_t gradient_center_opacity = 255;
uint8_t gradient_edge_opacity = 0;

size_t grid_size = 5;
size_t grid_width = grid_size*2;
size_t grid_height = grid_size*2;

float spacing = 2;

float boiler_width = 2*gradient_radius*grid_width*spacing;
float boiler_height = 2*gradient_radius*grid_height*spacing;

float window_width = boiler_width;
float window_height = boiler_height;

float MODEL_SCALE = 100.0;


float T_MAX =  573.15 * MODEL_SCALE;
float T_MIN =  273.15 * MODEL_SCALE;


size_t ATOM_NUMBER = grid_width*grid_height;


// 3*k/m
float KM = 0.3;

float g_acc = 5*9.81/MODEL_SCALE;

float DT = 1;

bool show_density = false;
size_t max_fps = 120;

bool debug_mode = false;


////////////////////////////////////////////////////////////////////////////////


std::string to_string(const Vector2& v)
{
    return "{ .x=" + std::to_string(v.x) + ", .y=" + std::to_string(v.y) + " }";
}

struct Atom
{
    Vector2 pos;
    Vector2 vel;
    Vector2 acc;
    float t = 0;
};

std::string to_string(const Atom& a)
{
    return "{ .pos=" + to_string(a.pos)
        + ", .vel=" + to_string(a.vel)
        + ", .acc=" + to_string(a.acc)
        + ", .t=" + std::to_string(a.t) + " }";
}


int sign(float x)
{
	return x >= 0 ? 1 : -1;
}

Vector2 t_to_vel(float t, Vector2 vel)
{
    float x2 = vel.x*vel.x;
    float y2 = vel.y*vel.y;
    float s = x2 + y2;
    float qx = s > 0 ? x2 / s : 0;
    float qy = s > 0 ? y2 / s : 0;
    return { sign(vel.x) * sqrt(t*KM*qx),
             sign(vel.y) * sqrt(t*KM*qy) };
}


float vel_to_t(Vector2 vel)
{
    return (vel.x*vel.x + vel.y*vel.y)/KM;
}


float model_width = boiler_width*MODEL_SCALE;
float model_height = boiler_height*MODEL_SCALE;

float model_radius = gradient_radius*MODEL_SCALE;


////////////////////////////////////////////////////////////////////////////////

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

auto grid = new_grid(grid_width, grid_height);

void print_grid(const grid_t& g)
{
    size_t iy = 0;
    for (const auto& row : g) {
        size_t ix = 0;
        for (const auto& cell : row) {
            if (cell.size() > 0) {
                printf("Cell %lu %lu:\n", iy, ix);
                for (const Atom* atom : cell) {
                    printf("    %p: %s\n", atom, to_string(*atom).c_str());
                }
            }
            ix++;
        }
        iy++;
    }
}

void index_atom(grid_t& g, Atom& atom)
{
    size_t ix = size_t(atom.pos.x / MODEL_SCALE / (2*gradient_radius*spacing));
    size_t iy = size_t(atom.pos.y / MODEL_SCALE / (2*gradient_radius*spacing));
    g[iy][ix].push_back(&atom);
}

std::vector<Atom> atoms;

void print_atoms()
{
    for (const auto& atom : atoms) {
        printf("%p: %s\n", &atom, to_string(atom).c_str());
    }
}

Atom& create_atom(Vector2 pos, Vector2 v, Vector2 a)
{
    atoms.emplace_back(Atom{pos, v, a, vel_to_t(v)});
    return atoms.back();
}

Atom& create_test_atom(float ix, float iy, float dx, float dy)
{
    float x = ix*2*gradient_radius*spacing + gradient_radius*spacing;
    float y = iy*2*gradient_radius*spacing + gradient_radius*spacing;
    x *= MODEL_SCALE;
    y *= MODEL_SCALE;
    return create_atom({x, y}, {dx, dy}, {0, 0});
}

void setup_debug_mode()
{
    // g_acc = 0
    show_density = true;
    max_fps = 120;

    atoms.reserve(4);

    index_atom(grid, create_test_atom(1.5, 0.5, -100, 100));
    index_atom(grid, create_test_atom(4, 2.5, 100, -100));
    //print_atoms();
    //print_grid(grid);
}

std::random_device rand_dev;
std::default_random_engine dre1(rand_dev());

int randrange(int from, int to)
{
    std::uniform_int_distribution<int> uniform_dist(from, to);
    return uniform_dist(dre1);
}

void setup()
{
    atoms.reserve(grid_width*grid_height);

    size_t iy = 0;
    for (auto& row : grid) {
        size_t ix = 0;
        for (auto& cell : row) {
            float x = ix*2*gradient_radius*spacing + gradient_radius*spacing;
            float y = iy*2*gradient_radius*spacing + gradient_radius*spacing;
            x *= MODEL_SCALE;
            y *= MODEL_SCALE;

            float t = randrange(int(T_MIN), int(T_MAX));
            float qx = randrange(-10, 10);
            float qy = randrange(-10, 10);

            Atom& atom = create_atom({x, y}, t_to_vel(t, {qx, qy}), {0, -g_acc});

            cell.push_back(&atom);

            ix++;
        }
        iy++;
    }
}


////////////////////////////////////////////////////////////////////////////////


Vector2 model_pos_to_graphics(Vector2 p)
{
    float x = p.x/MODEL_SCALE - window_width/2 - gradient_radius;
    float y = -p.y/MODEL_SCALE + window_height/2 - gradient_radius;
    return {x, y};
}

uint8_t temperature_to_sprite_index(float t)
{
    float T = T_MAX - T_MIN;
    float t_min = T_MIN - T;
    float t_max = T_MAX + T;
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


////////////////////////////////////////////////////////////////////////////////


float margin_h = (window_width - boiler_width)/2;
float margin_v = (window_height - boiler_height)/2;

Color gradient_inner_color(uint8_t i)
{
    return {uint8_t(i), 0, uint8_t(255-i), gradient_center_opacity};
}

Color gradient_outer_color(uint8_t i)
{
    return {uint8_t(i), 0, uint8_t(255-i), gradient_edge_opacity};
}

std::vector<Texture2D> atom_colors;

std::vector<Texture2D> generate_gradients()
{
    std::vector<Texture2D> textures;
    for (size_t i = 0; i < 256; i++) {
        Image radial_gradient = GenImageGradientRadial(gradient_radius * 2,
                                                       gradient_radius * 2,
                                                       0.0f,
                                                       gradient_inner_color(i),
                                                       gradient_outer_color(i));
        textures.push_back(LoadTextureFromImage(radial_gradient));
        UnloadImage(radial_gradient);
    }
    return textures;
}


////////////////////////////////////////////////////////////////////////////////


float distance2(Vector2 a, Vector2 b)
{
    return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}

float dot(Vector2 a, Vector2 b)
{
    return a.x*b.x + a.y*b.y;
}

Vector2 norm(Vector2 a)
{
    float s = sqrt(dot(a, a));
    if (s != 0)
        return {a.x/s, a.y/s};
    else
        return {0, 0};
}

Vector2 proj(Vector2 a, Vector2 b)
{
    float s = dot(a, b)/dot(b, b);
    return {b.x*s, b.y*s};
}

Vector2 vsum(Vector2 a, Vector2 b)
{
    return {a.x + b.x, a.y + b.y};
}

Vector2 vdif(Vector2 a, Vector2 b)
{
    return {a.x - b.x, a.y - b.y};
}

Vector2 vmul(Vector2 a, float k)
{
    return {a.x*k, a.y*k};
}

Vector2 invert(Vector2 a)
{
    return {-a.x, -a.y};
}

Vector2 orto(Vector2 a)
{
    return {-a.y, a.x};
}


////////////////////////////////////////////////////////////////////////////////

float model_d2 = 4*model_radius*model_radius;


void impact(Atom& atom1, Atom& atom2, float dt)
{
    float d2 = distance2(atom1.pos, atom2.pos);
    if (d2 > model_d2) {
        return;
    }

    // print("before:", atom1.vel.x, atom1.vel.y, atom2.vel.x, atom2.vel.y)
    // print("before:", atom1.acc.x, atom1.acc.y, atom2.acc.x, atom2.acc.y)

    // vector from center of atom1 to center of atom2 goes through impact place
    Vector2 s12 = vdif(atom2.pos, atom1.pos);
    Vector2 s21 = invert(s12);

    // print(s12.x, s12.y)

    // normal vectors to impact line - atoms will exchange velocity projections on this vector
    Vector2 n12 = norm(s12);
    Vector2 n21 = invert(n12);

    // print(n12.x, n12.y, n21.x, n21.y)

    // velocity of atom1 in projection on vector from atom1 to atom2 - then goes to atom2
    Vector2 v1n12 = proj(atom1.vel, s12);
    // velocity of atom2 in projection on vector from atom2 to atom1 - then goes to atom1
    Vector2 v2n21 = proj(atom2.vel, s21);
    // left velocity - saved
    Vector2 v1p = vdif(atom1.vel, v1n12);
    Vector2 v2p = vdif(atom2.vel, v2n21);

    // exchange impulse
    atom1.vel = vsum(v1p, v2n21);
    atom2.vel = vsum(v2p, v1n12);

    float d = sqrt(d2);
    float r = model_radius - d/2;
    r += 1.0;

    atom1.pos = vsum(atom1.pos, vmul(n21, r));
    atom2.pos = vsum(atom2.pos, vmul(n12, r));

    // print("after:", atom1.vel.x, atom1.vel.y, atom2.vel.x, atom2.vel.y)
    // print("after:", atom1.acc.x, atom1.acc.y, atom2.acc.x, atom2.acc.y)
}


grid_t update_model(const grid_t& g, float dt)
{
    grid_t result_grid = new_grid(grid_width, grid_height);

    int iy = 0;
    for (const auto& row : g) {
        int ix = 0;
        for (const auto& cell : row) {
            for (Atom* atom : cell) {

                for (int pix = ix-1; pix < ix+1; pix++) {
                    for (int piy = iy-1; piy < iy+1; piy++) {
                        if (pix >= 0 && piy >= 0 && pix < int(row.size()) && piy < int(grid.size())) {
                            for (Atom* a : grid[piy][pix]) {
                                if (a != atom) {
                                    impact(*atom, *a, dt);
                                }
                            }
                        }
                    }
                }

                atom->vel.x += atom->acc.x*dt;
                atom->vel.y += atom->acc.y*dt;
                atom->t = vel_to_t(atom->vel);

                atom->pos.x += atom->vel.x*dt;
                atom->pos.y += atom->vel.y*dt;

                // borders
                if (atom->pos.x < model_radius) {
                    atom->pos.x = model_radius+1;
                    atom->vel.x = -atom->vel.x;
                }
                if (atom->pos.x > model_width-model_radius) {
                    atom->pos.x = model_width-model_radius-1;
                    atom->vel.x = -atom->vel.x;
                }

                if (atom->pos.y < model_radius) {
                    atom->pos.y = model_radius+1;
                    atom->vel.y = -atom->vel.y;
                    atom->t = T_MAX;
                    atom->vel = t_to_vel(atom->t, atom->vel);
                }
                if (atom->pos.y > model_height-model_radius) {
                    atom->pos.y = model_height-model_radius-1;
                    atom->vel.y = -atom->vel.y;
                    atom->t = T_MIN;
                    atom->vel = t_to_vel(atom->t, atom->vel);
                }

                index_atom(result_grid, *atom);
            }
            ix++;
        }
        iy++;
    }

    return result_grid;
}


void render(const grid_t& g)
{
    if (show_density) {
        size_t max_atoms_in_cell=spacing*spacing*2;

        int iy = 0;
        for (const auto& row : g) {
            int ix = 0;
            for (const auto& cell : row) {

                int c = int(cell.size()*255/max_atoms_in_cell);
                uint8_t u = c <= 255 ? 255-c : 0;
                Color color = {u, u, u, 255};

                float r_w = 2*gradient_radius*spacing;
                float r_h = 2*gradient_radius*spacing;
                Rectangle r{
                    ix*r_w - window_width/2,
                    iy*r_h - window_width/2,
                    r_w,
                    r_h
                };
                DrawRectangleRec(r, color);

                ix++;
            }
            iy++;
        }

    } else {
        DrawRectangleRec({-window_width/2, -window_height/2, window_width, window_height}, LIGHTGRAY);
    }

    for (const auto& atom : atoms) {
        Texture2D& tx = atom_colors[temperature_to_sprite_index(atom.t)];
        Vector2 pos = model_pos_to_graphics(atom.pos);
        DrawTextureEx(tx, pos, 0, 1.0, WHITE);
    }

    //str = "%.1f fps, %.2f/%.2f" % (last_fps, last_model_time, last_render_time)
}


////////////////////////////////////////////////////////////////////////////////


int main(void)
{
    if (debug_mode) {
        setup_debug_mode();
    } else {
        setup();
    }

    const int screenWidth = window_width;
    const int screenHeight = window_height;

    InitWindow(screenWidth, screenHeight, "2D atomic gas");

    float border_thickness = 10;
    auto blocks = Box({-300, -150}, 600, 300, GRAY, border_thickness);

    // Draw 256 gradients before start
    atom_colors = generate_gradients();

    Camera2D camera = { {0} };
    camera.target = {0, 0};
    camera.offset = (Vector2){ screenWidth / 2.0f, screenHeight / 2.0f };
    camera.rotation = 0.0f;
    camera.zoom = 1.0f;

    SetTargetFPS(max_fps);

    while (!WindowShouldClose()) {

        grid = update_model(grid, DT);

        BeginDrawing();
        ClearBackground(GREEN);
        BeginMode2D(camera);

        //for (const auto& block : blocks)
        //    DrawRectangleRec(block.rect, block.color);

        render(grid);

        EndMode2D();
        DrawFPS(20, 5);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}
