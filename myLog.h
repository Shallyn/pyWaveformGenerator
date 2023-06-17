/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_MYLOG__
#define __INCLUDE_MYLOG__

#include <string.h>
 


/*通用字符串存储大小定义*/
#define STR_COMM_SIZE 128
#define STR_MAX_SIZE 1024

#define MAX_LOG_FILE_NUM		  (99)
 
#define NUMBER(type) sizeof(type)/sizeof(type[0])
 
#define __FILENAME__ (strrchr(__FILE__, '/') ? (strrchr(__FILE__, '/') + 1) : __FILE__)
 
/*日志类型*/
enum
{
    LOG_0 = 0,
    LOG_CRITICAL,/*调试日志*/
    LOG_ERROR,/*错误日志*/
    LOG_WARNING,/*告警日志*/
    LOG_INFO,/*运行日志*/
    LOG_DEBUG,/*系统日志*/
    BUTTOM
};

typedef struct tagCtrlParams
{
    int level;
    char flog[STR_COMM_SIZE];
    char fdebug[STR_COMM_SIZE];
}CtrlParams;


/*将日志输出到终端*/
#define PRINT_LOG_TO_TERM (0)
/*将日志输出到文件中*/
#define PRINT_LOG_TO_FILE (1)
 
/*调试日志宏定义*/
#define DEBUG_PRINT 0
#define LOG_PRINT(fmt, ...) do{\
	if(DEBUG_PRINT)\
	{\
		fprintf(stderr, fmt"  [line:%d] [%s]\n", ##__VA_ARGS__, __LINE__, __FUNCTION__);\
	}\
}while(0);
 
/*错误日志打印(在日志打印模块还未启动时使用)*/
#define LOG_ERR(fmt, ...) do{\
	fprintf(stderr, "[ERROR]  "fmt"  [line:%d] [%s]\n", ##__VA_ARGS__, __LINE__, __FUNCTION__);\
}while(0);
 
 
 
/*存储日志标记. 0-不存储日志, 1-存储日志*/
extern unsigned long g_ulPrintLogPlaceFlag;
 
/*根据等级判断是否打印调试日志标记*/
extern unsigned long g_ulPrintDebugLogFlag;
 
unsigned long LOG_PrintLog(unsigned char ucType, unsigned char *pucLogInfo);
 
 
/*日志打印宏定义*/
#define PRINT_LOG_INFO(type, fmt, ...) do{\
    if((g_ulPrintDebugLogFlag >= (long) type)) \
    {\
        unsigned char ucLogInfo[STR_MAX_SIZE] = {0}; \
        snprintf((char *)ucLogInfo, sizeof(ucLogInfo) - 1, fmt"  [%s] [line:%d] [%s]\n", ##__VA_ARGS__, __FILENAME__, __LINE__, __FUNCTION__); \
        LOG_PrintLog(type, ucLogInfo); \
    }\
}while(0)
 
 
/*是否打印调试日志标记,0-不打印调试日志,1-打印调试日志*/
extern void LOG_SetPrintDebugLogFlag(unsigned long flag);
/*存储日志标记. 0-不存储日志, 1-存储日志*/
extern void LOG_SetPrintLogPlaceFlag(unsigned long flag);
 
extern unsigned long LOG_Init(char* ucLogFileName, unsigned long ulFileSize);
extern void LOG_Destroy(void);

#endif

