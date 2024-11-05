
typedef struct{

    int N;      // 表序号
    int Ncol;   // 表列数
    int Nrow;   // 表行数
    double **values;    // 二维表内容

} TABLE;

typedef struct{

    int N;          // 表总数
    TABLE** tables; // 表指针

} TABLELIST;


char getWord(FILE* ff, char* s);

void tableMalloc(TABLE* t);

TABLE* tableRead(FILE* ff);

void tableDisplay(TABLE* t);

TABLELIST* fReadTables(char* fileName);

void readTablesFree(TABLELIST* tl);

/*

int main()
{

    TABLELIST* tl;

    tl = fReadTables("./mesh.csv");    

    tableDisplay(tl->tables[1]);

    return 0;

}

*/
