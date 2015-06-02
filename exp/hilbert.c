#include <stdio.h>

typedef unsigned int uint;


void rot(uint width, uint * x, uint * y, uint rx, uint ry)
{
    if (ry == 0U)
    {
        if (rx == 1U)
        {
            *x = width - 1U - *x;
            *y = width - 1U - *y;
        }
        uint t = *x;
        *x = *y;
        *y = t;
    }
}


uint xy_to_d (uint width, uint x, uint y)
{
    uint rx;
    uint ry;
    uint d = 0U;
    for (uint s = width / 2U; s > 0U; s /= 2U)
    {
        rx = (x & s) > 0U;
        ry = (y & s) > 0U;
        d += s * s * ((3U * rx) ^ ry);
        rot(s, &x, &y, rx, ry);
    }
    return d;
}


void d_to_xy (uint width, uint d, uint * x, uint * y)
{
    *x = *y = 0U;
    uint t = d;
    uint rx;
    uint ry;
    for (uint s = 1U; s < width; s *= 2U)
    {
        rx = 1U & (t / 2U);
        ry = 1U & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4U;
    }
}


int main(int argc, char *argv[])
{
    uint width = 13U;
    for (unsigned int i = 0U; i < width; i++)
    {
        for (unsigned int j = 0U; j < width; j++)
        {
            uint d = xy_to_d (width, i, j);
            printf ("%4u", d);
        }
        printf ("\n");
    }
    return 0;
}
