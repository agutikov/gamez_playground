#pragma once

#include <random>
#include <utility>
#include <vector>
#include <numeric>
#include <tuple>
#include <limits>
#include <fmt/core.h>
#include "raylib.h"


float randrange(float from, float to)
{
    static std::random_device rand_dev;
    static std::default_random_engine dre(rand_dev());

    std::uniform_real_distribution<float> uniform_dist(from, to);
    return uniform_dist(dre);
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


template<typename T>
struct stat_value
{
    struct stats_t
    {
        T min = 0;
        T mean = 0;
        T max = 0;

        operator std::string() const
        {
            return fmt::format("[{:.2f}, {:.2f}, {:.2f}]", min, mean, max);
        }
    };

    void consume(T v)
    {
        if (v < min) {
            min = v;
        }
        if (v > max) {
            max = v;
        }
        sum += v;
        count++;
    }
    stats_t produce()
    {
        if (count == 0) {
            return {};
        }
        stats_t s{min, sum/count, max};
        reset();
        return s;
    }

    stat_value()
    {
        reset();
    }

private:
    void reset()
    {
        min = std::numeric_limits<T>::max();
        max = std::numeric_limits<T>::min();
        sum = 0;
        count = 0;
    }

    T min;
    T max;
    T sum;
    size_t count;
};
