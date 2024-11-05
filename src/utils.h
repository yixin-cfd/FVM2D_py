
double** tableMallocDouble(int Nrow, int Ncol);  /*!< \brief 分配Nrow*Ncol 大小二维表格*/
 
int** tableMallocInt(int Nrow, int Ncol);

void tableFreeDouble(double** M, int Nrow);

void tableFreeInit(int** M, int Nrow);

double duration(struct timeval start, struct timeval stop);

