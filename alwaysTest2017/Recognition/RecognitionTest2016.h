#pragma once

//---------------------------------------------
// �㷨�ӿ�2
// DataBuf: 	[in] double buf[10240] 	һ���𶯵�����
// EnvironTemp: [in] double buf[820]	����ģ��
// IntrudeTemp: [in] double buf[820]	����ģ��
// Param:	[in] __para		4���ɵ��㷨����
//
// return: 0������ 1����
// �޸Ĳ���
// ���ǵ��Ժ����ݵ��������������2������
// int nDataSize	ÿ�����ݵĳ���
// __para &Param	���Լල��Ҫ�˽�Ĳ��ֲ���
//---------------------------------------------
typedef struct __PARA
{
	double fDistanceCtrl;	// ���Ҿ�����Сֵ
	double fReduceTail;		// <ReduceTail��Ϊ0
	int nJudgement_1;		// �����ж�1
	int nJudgement_2;		// �����ж�2
	// ��ӵ������չʾ����
	int ShowOutPos1;
	int ShowOutPos2;
	int DataEnd;	
}__para;
extern __para g_Param;

// �㷨����
//�ź�Ԥ����
double* DataFill(double *data_source, int nDataSize,int NFFT);	
//FFT
void DataFFT(double data[], int nn);		// ���Ż�
//FFT Optimization and Data Normalization
void DataCut(double X[], int nn, int Nx);
//ƽ���Ϳ�����
double sum_sqrt(double data[], int nn);	
//ģ��ƽ���Ϳ�����
double model_sqrt(double data[], int nn);	
//�������Ҿ���
double caculate_value(double data[], double data_sqrts,double models[], double model_sqrts, int ModelSize);
//�Ƚ��ź���ģ��֮������С
int max_data(double models[]);
//�źŷ���
int classify(double data[], double Environ[], double Intrude[], int ModelSize,int& ShowOut);
//FFT Reduce Tail and Data Normalization
void DataReduce(double X[], int nn, int Nx,double Limit);
// �ҳ�������ʼ��
int DataZeroJudge(double X[], int nn, int Nx);
// ����ж�ֵ
int GetInterval(int nn, int judge1, int judge2);
//�ж��ź�����
int choose_result(int chooses);

// �ӿ�2:   int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, __para &Param);
// ԭ�ӿ�:	int ProcessOneData(double* data_in,int nDataSize,int& nOut);
int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, int nDataSize);
int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, int nDataSize, int ModelSize, __para &Param);