#include <stdio.h>

void print_table(unsigned int width, unsigned int height)
{
    for (unsigned int i = 0; i < width; i++)
    {
        for (unsigned int j = 0; j < width; j++)
        {
            printf("%u", (i & 1U) ^ (j & 1U));
            unsigned int index = i * width + j;
            printf("%u ", ((index >> 3) & 1U) ^ (index & 1U));
        }
        printf("\n");
    }
}

void print_bits(unsigned int value)
{
    unsigned short size = 8U * sizeof(unsigned int);
    for (unsigned short i = 0; i < size; i++)
    {
        printf("%u", (value >> (size - i - 1)) & 1U);
    }
    printf(" :%u\n", value);
}


void deindex(unsigned int index, unsigned int width,
             unsigned int * x, unsigned int * y)
{
    *x = index % width;
    *y = index / width;
}


int main(int argc, char *argv[])
{
    print_table(10, 10);
    printf("\n\n");
    print_table(8, 8);
    printf("\n\n");
    print_table(11, 11);

    for (unsigned short i = 0; i < 21; i++)
    {
        printf("%u %u\n", i, 1U << i);
    }

    for (unsigned int i = 0; i < 1024; i++)
    {
        print_bits(1024U);
        print_bits(1024U * i);
    }

    unsigned short flag = 1U;

    for (unsigned int i = 0; i < 1024U * 1024U; i++)
    {
        unsigned int a = i % 1024U;
        unsigned int b = i & 1023;
        if (a != b)
        {
            flag = 0U;
        }
    }

    if (flag != 1U)
    {
        printf("*** EROR: You're thinking it wrong ***");
    }

    unsigned int size = 8U;
    unsigned int x;
    unsigned int y;

    for (unsigned int i = 0U; i < size * size; i++)
    {
        deindex(i, size, &x, &y);
        if (x == 0U) printf("\n");
        printf("%u ", (x == 0U) || ((x + 1) == size) || (y == 0U) || ((y + 1) == size));
    }

    return 0;
}
