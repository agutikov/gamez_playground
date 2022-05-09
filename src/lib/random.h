#pragma once

#include <random>
#include <utility>


float randrange(float from, float to)
{
    static std::random_device rand_dev;
    static std::default_random_engine dre1(rand_dev());

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
