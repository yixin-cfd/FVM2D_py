#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"

INPUT *inputInit(char *fileName, int N)
{
    INPUT *input = malloc(sizeof(INPUT));   // 给input分配内存
    FILE *ff = fopen(fileName, "r");        // 打开配置文件
    char c;                                 // 存储读取到的字符
    int name, value;                        // 指示当前读取的字符串是值还是命令名称
    int ii, jj;

    input->Nmax = N;                        // 最大字符串长度

    input->name = (char **)malloc(N * sizeof(char *)); // 为存储配置命令名称变量分配内存
    for (ii = 0; ii < N; ii++)
    {
        input->name[ii] = (char *)malloc(N * sizeof(char)); // 初始化
        input->name[ii][0] = '\0';
    }

    input->value = (char **)malloc(N * sizeof(char *)); // 为存储配置命令值变量分配内存
    for (ii = 0; ii < N; ii++)
    {
        input->value[ii] = (char *)malloc(N * sizeof(char)); // 初始化   
        input->value[ii][0] = '\0';
    }

    name = 1;           // 当前字符串是命令                 
    value = 0;          // 当前字符串不是值
    ii = 0;             // 第几个控制命令
    jj = 0;             // 字符串的待存字符索引
    c = fgetc(ff);      // 读取字符串
    while (c != EOF)    // 未到文件结尾
    {
        if ((c != ' ') & (c != '\r'))   // 有意义的字符
        {

            if (c == ',')               // 命令-值 分割符
            {
                name = 0;               // 当前字符串不是命令名称
                value = 1;              // 当前字符串是值
                input->name[ii][jj] = '\0'; // 给命令名加入'\0'
                jj = 0;                 // 新的值字符应该存放的索引
            }
            else if (c == '#')          // 注释语句
            {
                name = 0;
                value = 0;
                input->value[ii][jj] = '\0';
            }
            else if (c == '\n')         // 换行
            {
                if (value == 1)         // 值读取完毕
                {
                    name = 1;           // 读取新命令名称标记
                    value = 0;          // 清除读取值标记
                    input->value[ii][jj] = '\0';    // 上一次读取的值构成字符串
                    ii += 1;            // 新命令
                    jj = 0;             // 重置待存字符串位置索引
                }
                if (name == 1 && jj > 0)
                {
                    input->name[ii][jj] = '\0';
                    ii += 1;
                    jj = 0;
                }
            }
            else if (name)      // 命令名
            {
                input->name[ii][jj] = c;
                jj += 1;
            }
            else if (value)     // 值
            {
                input->value[ii][jj] = c;
                jj += 1;
            }
        }

        c = fgetc(ff);
    }

    fclose(ff);

    input->N = ii;      // 命令数量

    return input;
}

void inputPrint(INPUT *input)
{
    for (int ii = 0; ii < input->N; ii++)
    {
        printf(" %s, %s\n", input->name[ii], input->value[ii]);
    }
}

int inputNameIsInput(INPUT *input, char *name)
{

    int found = 0;

    for (int ii = 0; ii < input->N; ii++)
    {
        if (strcmp(name, input->name[ii]) == 0)
        {
            found = 1;
        }
    }

    return found;
}

char *inputGetValue(INPUT *input, char *name)
{
    char *s;
    int found = 0;

    for (int ii = 0; ii < input->N; ii++)
    {
        if (strcmp(name, input->name[ii]) == 0)
        {
            s = input->value[ii];
            found = 1;
        }
    }

    if (found == 0)
    {
        printf("Error: Value not found in the input file: %s.\n", name);
        exit(0);
    }

    return s;
}

void inputFree(INPUT *input)
{

    for (int ii = 0; ii < input->N; ii++)
    {
        free(input->name[ii]);
        free(input->value[ii]);
    }
    free(input->name);
    free(input->value);
}

/*

int main()
{

    INPUT* input = initInput("./input.csv");

    inputPrint(input);

    printf("\n %s \n", inputGetValue(input, "p"));

    return 0;

}

*/
