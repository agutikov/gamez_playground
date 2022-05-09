#pragma once

#include "raymath.h"
#include <string>
#include <sstream>
#include <iomanip>


Vector2 operator-(const Vector2& a)
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

