/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <pthread.h>
#include "myLog.h"


/*存储日志的文件名*/
 
static unsigned char g_ucLogFileName[MAX_LOG_FILE_NUM][STR_COMM_SIZE] = {{0}};
 
/*指明是g_ucLogFileName中的哪个文件*/
static unsigned char g_ucLogFileNo = 0;
 
/*输出日志位置标记,0-输出到终端,1-输出到日志文件*/
unsigned long g_ulPrintLogPlaceFlag = 0;
/*是否打印调试日志标记,0-不打印调试日志,1-打印调试日志*/
unsigned long g_ulPrintDebugLogFlag = 0;

/*日志文件大小*/
static unsigned long g_ulLogFileSize = 0;
 
/*日志文件句柄*/
static FILE* pFile = NULL;
 
/*日志存储互斥锁*/
static pthread_mutex_t g_stSaveLogMutexLock;
 
/*日志模块初始化标记*/
static unsigned long g_ulLogInitFlag = 0;
 
 
void LOG_SetPrintLogPlaceFlag(unsigned long flag)
{
	g_ulPrintLogPlaceFlag = flag;
}
 
void LOG_SetPrintDebugLogFlag(unsigned long flag)
{
	g_ulPrintDebugLogFlag = flag;
}
 
/*****************************************************************
** 函数名: get_file_size
** 输　入: char *path
** 输　出:
** 功能描述:获取指令文件大小
** 返回值: long
****************************************************************/
static long get_file_size(const char *path)
{
	long filesize = -1;
	struct stat statbuff;
	
	if(stat(path, &statbuff) < 0){
	    return filesize;
	}
	else{
	    filesize = statbuff.st_size;
	}
	return filesize;
}
 
/*****************************************************************
** 函数名: unsigned long LOG_PrintLogTime
** 输　入:  unsigned long ulBufLen 存储时间的空间长度
** 输　出:unsigned char *ucTime  存储时间
** 功能描述:日志输出
** 返回值:unsigned long
****************************************************************/
unsigned long LOG_PrintLogTime(unsigned char *ucTime, unsigned long ulBufLen)
{
    struct tm* pstTmSec;
    struct timeval stTmMsec;
	
	if(NULL == ucTime)
	{
		return -1;
	}
	gettimeofday(&stTmMsec, NULL);
	pstTmSec = localtime(&stTmMsec.tv_sec);
	snprintf((char *)ucTime, ulBufLen - 1, "%04d-%02d-%02d %02d:%02d:%02d %03dms",
            pstTmSec->tm_year + 1900, pstTmSec->tm_mon + 1, pstTmSec->tm_mday, pstTmSec->tm_hour, 
            pstTmSec->tm_min, pstTmSec->tm_sec, stTmMsec.tv_usec / 1000);
	
	return 0;
}
 
/*****************************************************************
** 函数名: 	unsigned char *LOG_LogTypeToStr
** 输　入:  	unsigned char ucType  日志类型
			unsigned long ulBufLen 存储日志类型字符串空间的长度
** 输　出:unsigned char *pucTypeString 根据日志类型将其转换成相应的字符串
** 功能描述:根据日志类型转换成相应的字符串
** 返回值:unsigned long
****************************************************************/
unsigned long LOG_LogTypeToStr(unsigned char ucType, unsigned char *pucTypeString, unsigned long ulBufLen)
{
	if(NULL == pucTypeString)
	{
		return -1;
	}
	/*防止发生越界*/
	ulBufLen -= 1;
 

    if(PRINT_LOG_TO_TERM == g_ulPrintLogPlaceFlag)
    {
        switch(ucType)
        {
            case LOG_DEBUG:
            {
                strncpy((char *)pucTypeString, "\033[1;32mDEBUG\033[0m", ulBufLen);
                break;
            }
            case LOG_ERROR:
            {
                strncpy((char *)pucTypeString, "\033[1;31mERROR\033[0m", ulBufLen);
                break;
            }
            case LOG_WARNING:
            {
                strncpy((char *)pucTypeString, "\033[1;33mWARNING\033[0m", ulBufLen);
                break;
            }
            case LOG_INFO:
            {
                strncpy((char *)pucTypeString, "\033[1;37mINFO\033[0m", ulBufLen);
                break;
            }
            case LOG_CRITICAL:
            {
                strncpy((char *)pucTypeString, "\033[47;31mCRITICAL\033[0m", ulBufLen);
                break;
            }
            default:
            {
                strncpy((char *)pucTypeString, "UNKNOWN", ulBufLen);
                break;
            }
        }

    }
    else
    {
        switch(ucType)
        {
            case LOG_DEBUG:
            {
                strncpy((char *)pucTypeString, "DEBUG", ulBufLen);
                break;
            }
            case LOG_ERROR:
            {
                strncpy((char *)pucTypeString, "ERROR", ulBufLen);
                break;
            }
            case LOG_WARNING:
            {
                strncpy((char *)pucTypeString, "WARNING", ulBufLen);
                break;
            }
            case LOG_INFO:
            {
                strncpy((char *)pucTypeString, "INFO", ulBufLen);
                break;
            }
            case LOG_CRITICAL:
            {
                strncpy((char *)pucTypeString, "CRITICAL", ulBufLen);
                break;
            }
            default:
            {
                strncpy((char *)pucTypeString, "UNKNOWN", ulBufLen);
                break;
            }
        }
    }
	return 0;
}
 
/*****************************************************************
** 函数名: unsigned long LOG_OpenLogFile
** 输　入:  void
** 输　出:void
** 功能描述:打开日志文件
** 返回值:unsigned long
****************************************************************/
unsigned long LOG_OpenLogFile(void)
{
	char *path = (char*)g_ucLogFileName[g_ucLogFileNo];
	char *flag = NULL;
	int len = 0;
	
	/*判断文件是否已经打开*/
	if(NULL != pFile)
	{
		LOG_PRINT("[ACTION] file opened!");
		return 0;
	}
	/*判断文件名是否有定义*/
	if(NULL == path)
	{
		LOG_PRINT("[ERROR] file name is NULL.");
		return -1;
	}
	
	/*判断文件是否存在*/
	if (!access(path, 0))
	{
		/*获取文件大小*/
		if (0 > (len = get_file_size(path)))
		{
			LOG_PRINT("[ERROR] get file size failed!");
			return -1;
		}
	}
	flag = (len > 0 && len < g_ulLogFileSize) ? "a" : "w";
	
	/*打开文件*/
	pFile = fopen(path, flag);
	if(NULL == pFile)
	{
		LOG_PRINT("[ERROR] open file failed!");
		return -1;
	}
	LOG_PRINT("[DEBUG] open file name = %s", path);
	return 0;
}
 
/*****************************************************************
** 函数名: LOG_PrintLog
** 输　入:  unsigned char *ucLogInfo  需要打印或者存储的日志信息
			unsigned char ucType 日志类型
** 输　出:void
** 功能描述:日志输出
** 返回值:unsigned long
****************************************************************/
unsigned long LOG_PrintLog(unsigned char ucType, unsigned char *pucLogInfo)
{
	unsigned long ulResult = 0;
	unsigned long ulFileLen = 0;
	unsigned char ucTime[STR_COMM_SIZE] = {0};
	unsigned char ucLogTypeStr[STR_COMM_SIZE] = {0};
	unsigned char ucLogInfo[STR_MAX_SIZE] = {0};
	
	if(NULL == pucLogInfo)
	{
		return -1;
	}
	
	/*将日志类型转换成字符串*/
	ulResult = LOG_LogTypeToStr(ucType, ucLogTypeStr, sizeof(ucLogTypeStr));
	/*获取生成日志的时间*/
	ulResult += LOG_PrintLogTime(ucTime, sizeof(ucTime));
	if(0 != ulResult)
	{
		return -1;
	}
	snprintf((char *)ucLogInfo, sizeof(ucLogInfo) - 1, "[%s] [%s] %s", ucTime, ucLogTypeStr, pucLogInfo);
	/*判断是否打印调试日志*/
	if(PRINT_LOG_TO_TERM == g_ulPrintLogPlaceFlag)
	{
		fprintf(stderr, "%s", ucLogInfo);
		return 0;
	}
	/*加锁保护文件操作*/
	pthread_mutex_lock(&g_stSaveLogMutexLock);
	/*打开日志文件*/
	(void)LOG_OpenLogFile();
	if(NULL != pFile)
	{
		fputs((char *)ucLogInfo, pFile);
		ulFileLen = ftell(pFile);
		LOG_PRINT("file len = %ld", ulFileLen);
		if(ulFileLen >= g_ulLogFileSize)
		{
			fclose(pFile);
			pFile = NULL;
			g_ucLogFileNo = (g_ucLogFileNo + 1) % MAX_LOG_FILE_NUM;
		}
	}
	pthread_mutex_unlock(&g_stSaveLogMutexLock);
	return 0;
}
 
/*****************************************************************
** 函数名: LOG_Init
** 输　入:  const unsigned char* ucLogFileName  用来保存日志的文件名
			unsigned long ulFileSize 存储日志的文件大小
** 输　出:void
** 功能描述:日志打印初始化
** 返回值:unsigned long
****************************************************************/
unsigned long LOG_Init(char* ucLogFileName, unsigned long ulFileSize)
{
	unsigned int i = 0;
	/*判断参数的合法性*/
	if((NULL == ucLogFileName) || !(ulFileSize > 0))
	{
		return -1;
	}
	/*判断是否将日志输出到日志文件*/
	if((PRINT_LOG_TO_FILE != g_ulPrintLogPlaceFlag) || (0 != g_ulLogInitFlag))
	{
		//fprintf(stderr, "g_ulPrintLogPlaceFlag = %ld g_ulLogInitFlag = %ld\n",g_ulPrintLogPlaceFlag,g_ulLogInitFlag);
		LOG_PRINT("print log to termination!!");
		return 0;
	}
	
	/*记录日志模块已经被初始化(防止改模块被重复初始化)*/
	g_ulLogInitFlag = 1;
	
	/*生成存储日志的文件名*/
	for(i = 0; i < NUMBER(g_ucLogFileName); i++)
	{
		snprintf((char *)g_ucLogFileName[i], sizeof(g_ucLogFileName[i]) - 1, "%s_%02d", ucLogFileName, i);
		LOG_PRINT("Log File: %s", g_ucLogFileName[i]);
		//fprintf(stderr, "Log File: %s\n", g_ucLogFileName[i]);
	}
	/*设置日志文件大小*/
	g_ulLogFileSize = ulFileSize;
	pthread_mutex_init(&g_stSaveLogMutexLock, NULL);
 
	return 0;
}
 
/*****************************************************************
** 函数名: LOG_Destroy
** 输　入:  void
** 输　出:void
** 功能描述:日志打印资源释放
** 返回值:unsigned long
****************************************************************/
void LOG_Destroy(void)
{
	if(pFile != NULL)
	{
		fclose(pFile);
		pFile = NULL;
	}
	pthread_mutex_destroy(&g_stSaveLogMutexLock);
	return;
}
