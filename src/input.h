
typedef struct
{

    int N;                  /*!< \brief 配置变量总数  */     
    int Nmax;               /*!< \brief 字符串最大长度  */   
    char** name;            /*!< \brief 存储操作名  */  
    char** value;           /*!< \brief 存储操作值  */  

} INPUT;

INPUT* inputInit(char* fileName, int N);        //!<  配置文件初始化

void inputPrint(INPUT* input);                  //!< 打印配置文件

int inputNameIsInput(INPUT* input, char* name); //!< 检查命令是否在配置文件中

char* inputGetValue(INPUT* input, char* name);  //!< 读取命令对应的值

void inputFree(INPUT* input);   //!< 析构函数
