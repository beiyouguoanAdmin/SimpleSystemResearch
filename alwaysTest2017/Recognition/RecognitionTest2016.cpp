#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <YangTest.h>

/**********************************************************************
���¶���ɵ�������
#define ReduceTail 0.08		// <0.08��Ϊ0���䶯��Χ��0.001~0.2
#define CosDistance 0.6		// ���Ҿ�����Сֵ,�䶯��Χ0.5~1
#define Judgement_1 1000	// �����ж�1,�䶯��ΧΪ100~1000
#define Judgement_2 2000	// �����ж�2,�䶯��ΧΪ500~3000
�̶���ZeroPoint
ɾ���ˣ�Judgement_3 ,StrictParameter_DOWN��StrictParameter_UP


��������Ҫ���������
ShowOutPos1չʾ����ڶ����㷨����������δ����Ϊ-1
ShowOutPos1����GetInterval()��max_data()���ۺ����
10���ڶ����㷨���ֶȲ���
100~200���ڶ����㷨����
200~300���ڶ����㷨������

ShowOutPos2չʾ��һ���㷨���ж��������������˵ڶ����㷨���-1
1000����һ���㷨������
100����һ���㷨����

DataEnd�����ֹλ��
һ������Ȼ������500����
��Ϊ������1500~3000����

pengkun��������2016.09.17
***********************************************************************/

#define PI	M_PI
#define NUM_OF_MODEL 2		//����ģ������
#define nZeroPoint 40		//�����������Ŀ
//�ź�����
//�ɿز�������
/*
// �����жϿɵ�����
#define ReduceTail 0.08		// <0.08��Ϊ0
#define Judgement_1 1000	// �����ж�1
#define Judgement_2 2000	// �����ж�2
// kmeanģ��ɵ�����
#define CosDistance 0.6		// ���Ҿ�����Сֵ
*/
__para g_Param; 

//�ź�Ԥ����
double* DataFill(double *data_source, int nDataSize, int NFFT)
{
	int i;
	//printf("NFFT = %d\n", NFFT);

	double* temp = new double[2 * NFFT + 1];

	/* Storing x(n) in a complex array to make it work with four1.
	This is needed even though x(n) is purely real in this case. */
	for (i = 0; i<nDataSize; i++)
	{
		temp[2 * i + 1] = data_source[i];
		temp[2 * i + 2] = 0.0;
	}
	/* pad the remainder of the array with zeros (0 + 0 j) */
	for (i = nDataSize; i<NFFT; i++)
	{
		temp[2 * i + 1] = 0.0;
		temp[2 * i + 2] = 0.0;
	}
	return temp;
}

//FFT
void DataFFT(double data[], int nn)
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	// reverse-binary reindexing
	n = nn << 1;	// nn=16384,2^14; n=2^15
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			tempr = data[j];     data[j] = data[i];     data[i] = tempr;
			tempr = data[j + 1]; data[j + 1] = data[i + 1]; data[i + 1] = tempr;
		}
		m = n >> 1;		// m = nn; // m=2^14
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	// here begins the Danielson-Lanczos section
	mmax = 2;
	while (n > mmax) {
		istep = 2 * mmax;
		theta = -(2.0*PI) / mmax;		
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr*data[j] - wi*data[j + 1];
				tempi = wr*data[j + 1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j + 1] = data[i + 1] - tempi;
				data[i] += tempr;
				data[i + 1] += tempi;
			}
			wtemp = wr;
			wr = wr*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}
}

//FFT Reduce Tail and Data Normalization
void DataReduce(double X[], int nn, int Nx,double Limit)
{
	int i; double tempMax = 0.0; double tempMin = 1000.0;
	for (i = 0; i<Nx; i++)
	{
		X[2 * i + 1] = 2 * (sqrt(X[2 * i + 1] * X[2 * i + 1] + X[2 * i + 2] * X[2 * i + 2])) / Nx;
		if (tempMax<X[2 * i + 1])
			tempMax = X[2 * i + 1];
		if (tempMin>X[2 * i + 1])
			tempMin = X[2 * i + 1];
	}
	for (i = 0; i<nn; i++)
	{
		X[2 * i + 1] = (X[2 * i + 1] - tempMin) / (tempMax - tempMin);
		if(X[2 * i + 1]< Limit)
			X[2 * i + 1]=0;
		//printf("X[%d] = %f\n", i, X[2 * i + 1]);
	}
}

// �ҳ�������ʼ��
int DataZeroJudge(double X[], int nn, int Nx)
{
	int xi = 0;
	int addi = 0;
	while (xi<Nx)
	{
		if (X[2 * xi + 1] == 0)
			addi++;
		else
			addi = 0;
		if (addi == nn)
			break;
		else
			//std::cout << addi << std::endl;
			xi++;
	}
	//std::cout << xi<<std::endl;
	return xi;
}

// ����ж�ֵ
int GetInterval(int nn, int judge1, int judge2)
{
	int location;
	if (nn < judge1)
		location = 1000;
	else if (nn <= judge2)
		location = 2;
	else
		location = 100;		// ###location = 1;
	return location;
}

//FFT Optimization and Data Normalization
void DataCut(double X[], int nn, int Nx)
{
	int i; double tempMax = 0.0; double tempMin = 1000.0;
	for (i = 0; i<Nx; i++)
	{
		X[2 * i + 1] = 2 * (sqrt(X[2 * i + 1] * X[2 * i + 1] + X[2 * i + 2] * X[2 * i + 2])) / Nx;
		if (tempMax<X[2 * i + 1])
			tempMax = X[2 * i + 1];
		if (tempMin>X[2 * i + 1])
			tempMin = X[2 * i + 1];
	}
	for (i = 0; i<nn; i++)
	{
		X[2 * i + 1] = (X[2 * i + 1] - tempMin) / (tempMax - tempMin);
		//printf("X[%d] = %f\n", i, X[2 * i + 1]);
	}
}
//ʵʱ����ƽ���Ϳ�����
double sum_sqrt(double data[], int nn)
{
	double data_r=0.0;
	for (int i = 1; i<nn; i++)		//��ȷ��i����ʼ
	{
		data_r += pow(data[2 * i + 1], 2.0);
	//	std::cout << data_r<<", "<<data[i]<<". ";
	}
	data_r = sqrt(data_r);
	return data_r;
}
//ģ��ƽ���Ϳ�����
double model_sqrt(double data[], int nn)
{
	double data_r=0.0;
	for (int i = 1; i<nn; i++)		//��ȷ��i����ʼ
	{
		data_r += pow(data[i], 2);
	//	std::cout << data_r<<", "<<data[i]<<". ";
	}
	data_r = sqrt(data_r);
	return data_r;
}

//�������Ҿ���
double caculate_value(double data[], double data_sqrts,double models[], double model_sqrts, int ModelSize)
{
	double add_result = 0.0;
	double cosdis = 0;
	int i;
	for (i = 1; i<ModelSize; i++)		//��ȷ��i����ʼ
	{
		add_result += data[2 * i + 1] * models[i];
	}
	cosdis = add_result/(data_sqrts * model_sqrts);    // 128*????????
	//printf("%d,%f,%d\t", add_result, data_sqrts, cosdis);
	return(cosdis);

}
//�Ƚ��ź���ģ��֮������С
int max_data(double models[NUM_OF_MODEL])
{
	int result = 10;	//### result = 0;
	int temp = 1;
	double max_of_data = models[0];
	if(models[0]<models[1])
	{
		max_of_data = models[1];
		temp=2;
	}
	if(max_of_data > g_Param.fDistanceCtrl)
	{
		result = 100* temp + 100 * max_of_data;
	}	
	return result;
}

//�ж��ź�����
// adjuestment
int choose_result(int chooses)
{
	int value = 0;
	if ( chooses > 200)	value = 10000;
	else if (chooses > 100)	value = 100;
	else value = 10;
	return value;
}

//�źŷ���
int classify(double data[], double Environ[], double Intrude[], int ModelSize, int& ShowOut)
{	
	
	int cos_result = 0;
	double value_cos[NUM_OF_MODEL];

	double data_sqrts = sum_sqrt(data, ModelSize);
	double Environ_sqrts = model_sqrt(Environ, ModelSize);
	double Intrude_sqrts = model_sqrt(Intrude, ModelSize);
	
	value_cos[0] = caculate_value(data, data_sqrts, Intrude, Intrude_sqrts, ModelSize);
	value_cos[1] = caculate_value(data, data_sqrts, Environ, Environ_sqrts, ModelSize);

	cos_result = max_data(value_cos);
	ShowOut = cos_result;
	cos_result = choose_result(cos_result);

	return cos_result;

}

int IntrudeProcess(double *DataBuf, double *EnvironTemp, double *IntrudeTemp, int nDataSize,int ModelSize, __para &Param)
{
	g_Param.fDistanceCtrl = Param.fDistanceCtrl;
	g_Param.fReduceTail = Param.fReduceTail;
	g_Param.nJudgement_1= Param.nJudgement_1;
	g_Param.nJudgement_2= Param.nJudgement_2;
	g_Param.ShowOutPos1 = -1;
	g_Param.ShowOutPos2 = -1;
	g_Param.DataEnd = -1;	
	
	int out; int NFFT;int data_end;int preout;	
	// �������ڳ����ж�
	double* data_pro1;			//�����ع�10240�����ݣ������¼������ж�
	// ��������kmean
	double* data_pro2;			//�����ع�10240�����ݣ�ʹ�����㷨ƥ��
	//double DataBuf[N*20];	//�����10240������
	//�ź�Ԥ��?
	//��ȡ����ʱ��FFT�ķ�Χ
	int ShowOut;
	NFFT = (int)pow(2.0, ceil(log((double)nDataSize) / log(2.0)));	
	data_pro1 = DataFill(DataBuf, nDataSize, NFFT);	
	DataFFT(data_pro1, NFFT);
	DataReduce(data_pro1, NFFT/2, NFFT/2, g_Param.fReduceTail);//ReduceTail);
	data_end=DataZeroJudge(data_pro1, nZeroPoint, NFFT/2 );//ZeroPoint, NFFT/2 );
	g_Param.DataEnd = data_end;
	delete[] data_pro1;
	preout=GetInterval(data_end, g_Param.nJudgement_1, g_Param.nJudgement_2);	//����������:100,2,1000

	Param.DataEnd = g_Param.DataEnd;	// ���չʾ
	
	if(preout == 2)
	{
		// ����kmean����
		data_pro2 = DataFill(DataBuf, nDataSize, NFFT);
		DataFFT(data_pro2, NFFT);
		DataCut(data_pro2, ModelSize, NFFT/2);	//֮ǰ�����������󣬵������ݴ�������Ӱ��
		//printf("\n\n");
		//�������Ҿ���
		out = classify(data_pro2, EnvironTemp, IntrudeTemp, ModelSize, ShowOut);
		g_Param.ShowOutPos1 = ShowOut;
		//�ж��ź�����
		delete[] data_pro2;

		Param.ShowOutPos1 = g_Param.ShowOutPos1;
		Param.ShowOutPos2 = g_Param.ShowOutPos2;
		
		if(ShowOut>100 && ShowOut<200) return 1;
		else return 0;
	}
	else
	{
		ShowOut=preout;
		g_Param.ShowOutPos2 = ShowOut;
		Param.ShowOutPos1 = g_Param.ShowOutPos1;
		Param.ShowOutPos2 = g_Param.ShowOutPos2;		
		
		if(preout==100) return 1;
		else return 0;
	//system("pause");
	}

	return 0;
}

