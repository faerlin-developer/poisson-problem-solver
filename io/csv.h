#pragma once

class Csv {
public:
    std::string name;

    explicit Csv(std::string name);

    void write(double *grid, int W, int H);
};
