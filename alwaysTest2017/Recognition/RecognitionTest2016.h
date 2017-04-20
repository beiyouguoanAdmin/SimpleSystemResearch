#pragma once

//---------------------------------------------
// 算法接口2
// DataBuf: 	[in] double buf[10240] 	一次震动的数据
// EnvironTemp: [in] double buf[820]	环境模板
// IntrudeTemp: [in] double buf[820]	入侵模板
// Param:	[in] __para		4个可调算法参数
//
// return: 0非入侵 1入侵
// 修改部分
// 考虑到以后数据的情况，所以增加2个参数
// int nDataSize	每组数据的长度
// __para &Param	测试监督需要了解的部分参数
//---------------------------------------------
typedef struct __PARA
{
	double fDistanceCtrl;	// 余弦距离最小值
	double fReduceTail;		// <ReduceTail归为0
	int nJudgement_1;		// 区域判定1
	int nJudgement_2;		// 区域判定2
	// 添加调配过程展示内容
	int ShowOutPos1;
	int ShowOutPos2;
	int DataEnd;	
}__para;
extern __para g_Param;

// 算法处理
//信号预处理
double* DataFill(double *data_source, int nDataSize,int NFFT);	
//FFT
void DataFFT(double data[], int nn);		// 可优化
//FFT Optimization and Data Normalization
void DataCut(double X[], int nn, int Nx);
//平方和开根号
double sum_sqrt(double data[], int nn);	
//模板平方和开根号
double model_sqrt(double data[], int nn);	
//计算余弦距离
double caculate_value(double data[], double data_sqrts,double models[], double model_sqrts, int ModelSize);
//比较信号与模板之间距离大小
int max_data(double models[]);
//信号分类
int classify(double data[], double Environ[], double Intrude[], int ModelSize,int& ShowOut);
//FFT Reduce Tail and Data Normalization
void DataReduce(double X[], int nn, int Nx,double Limit);
// 找出连零起始点
int DataZeroJudge(double X[], int nn, int Nx);
// 输出判定值
int GetInterval(int nn, int judge1, int judge2);
//判断信号类型
int choose_result(int chooses);

// 接口2:   int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, __para &Param);
// 原接口:	int ProcessOneData(double* data_in,int nDataSize,int& nOut);
int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, int nDataSize);
int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, int nDataSize, int ModelSize, __para &Param);