
typedef struct ELEM
{

    int ii;     /*!< \brief 单元索引 */ 
    int Np;     /*!< \brief 节点数量 */    
    int neiN;   /*!< \brief 邻接单元总数 */
    
    int* p;     /*!< \brief 存储节点索引 */
    int* f;
    
    double* P; 
    double d;       
    double omega;

    struct ELEM** neiL;    /*!< \brief 存储邻接单元 */
    
} ELEMENT;

/*!< \brief 边界段结构 */
typedef struct
{

    int Nelem;          /*!< \brief 边界单元数量 */
    
    char name[50];      /*!< \brief 边界段名称 */

    int flagBC;         /*!< \brief 边界条件类型 */

    ELEMENT** elemL;

} MESHBC;   

/*!< \brief 网格数据结构 */
typedef struct
{

    int Ndim;       /*!< \brief 维度 */
    int Nelem;      /*!< \brief 单元数量 */
    int Np;         /*!< \brief 节点数量 */        
    int Nmark;      /*!< \brief mark 数量 */   
    int Ncon;       /*!< \brief 共享边的数量 */
    int axi;
    int order;
    
    int** con;          /*!< \brief 共享边数组 */   

    double** p;         /*!< \brief 节点坐标数组 */

    ELEMENT** elemL;    /*!< \brief 单元数组 */
    
    MESHBC** bc;        /*!< \brief 存储边界段 */

} MESH;                 

typedef struct CON
{

    int data[4];
    struct CON* next; 

} CONNECTION;

/*!< \brief 边数据结构(链表节点) */
typedef struct FACE
{
    int full;   /*!< \brief 指示边是否已经记录过邻接单元 */
    int p0;     /*!< \brief 边的第一个节点 */
    int p1;     /*!< \brief 边的第二个节点 */
    int e0;     /*!< \brief 边所属单元 */
    int e1;     /*!< \brief 边邻接单元 */
    struct FACE* next;  /*!< \brief 下一条边 */

} FACETYPE;

/*!< \brief 读取关键字 */
char meshGetWord(FILE* ff, char* s);       
/*!< \brief 读取边界段 */
MESHBC* meshBCread(FILE* ff, int Nvar);     
/*!< \brief 网格初始化 */
MESH* meshInit(char* fileName, int Nvar, int axi);
/*!< \brief 打印边界信息 */
void meshPrintBC(MESHBC* bc);

void meshPrint(MESH* mesh);

void meshElementFree(ELEMENT* e);

void meshBCFree(MESHBC* bc);

void meshFree(MESH* mesh);

void elementCenter(ELEMENT* E, MESH* mesh, double* x, double* y);

double meshCalcOmegaTri(MESH* mesh, int p0, int p1, int p2);

double meshCalcDSlateral(MESH* mesh, int ii);

double meshCalcOmega(MESH* mesh, int ii);

void meshUpdateOmega(MESH* mesh);

double elementIsConnected(ELEMENT* e0, ELEMENT* e1, int* p0, int* p1);
/*!< \brief 用于判断两个边是否相同 */
int meshSameFace(FACETYPE* f0, FACETYPE* f1);

void meshCalcConnection1(MESH* mesh);
/*!< \brief 记录共享边信息 */
void meshCalcConnection2(MESH* mesh);

void meshPrintConnection(MESH* mesh, int N);

void meshCalcDS(MESH* mesh, int p0, int p1, double* dSx, double* dSy);

void meshCalcDS2(MESH* mesh, int p0, int p1, double* nx, double* ny, double* dS);

int meshBCIsConnect(ELEMENT* BCe, ELEMENT* e);
/*!< \brief 利用共享边信息更新单元数组的相应邻接单元 */
void meshCalcNeighbors(MESH* mesh);

void meshBCneighbors(MESHBC* bc, MESH* mesh);

double meshEdgeLength(MESH* mesh, int p0, int p1);

double meshMinEdge(MESH* mesh);

void meshPrintDStotal(MESH* mesh);

void meshCheckUse(MESH* mesh);

int meshPOri(MESH* mesh, ELEMENT* e, int p0, int p1);

void meshCheckBorderOrientation(MESH* mesh);

void meshFixBorderOrientation(MESH* mesh);

void meshCalcFaces(MESH* mesh);
