#pragma once

#include "raylib.h"
#include "raymath.h"
#include <vector>
#include <sstream>
#include <iomanip>


struct IntRectangle
{
    int x;
    int y;
    int width;
    int height;

    operator Rectangle() const
    {
        return {
            float(x),
            float(y),
            float(width),
            float(height)
        };
    }
};


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

std::ostream& operator<< (std::ostream& os, const Rectangle& r)
{
    os << std::fixed << std::setprecision(2)
        << "{.x=" << r.x << ", .y=" << r.y
        << " .width=" << r.width << ", .height=" << r.height
        << "}";
    return os;
}

std::string to_string(const Rectangle& r)
{
    std::stringstream ss;
    ss << r;
    return ss.str();
}

