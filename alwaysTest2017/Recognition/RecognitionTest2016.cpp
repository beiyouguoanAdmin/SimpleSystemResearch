#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <YangTest.h>

/**********************************************************************
重新定义可调参数：
#define ReduceTail 0.08		// <0.08归为0，变动范围是0.001~0.2
#define CosDistance 0.6		// 余弦距离最小值,变动范围0.5~1
#define Judgement_1 1000	// 区域判定1,变动范围为100~1000
#define Judgement_2 2000	// 区域判定2,变动范围为500~3000
固定了ZeroPoint
删除了：Judgement_3 ,StrictParameter_DOWN和StrictParameter_UP


建立我需要的输出变量
ShowOutPos1展示进入第二重算法的输出情况，未进入为-1
ShowOutPos1代表GetInterval()和max_data()的综合输出
10：第二重算法区分度不够
100~200：第二重算法报警
200~300：第二重算法不报警

ShowOutPos2展示第一重算法的判定情况，如果进入了第二重算法输出-1
1000：第一重算法不报警
100：第一重算法报警

DataEnd代表截止位置
一般让自然环境在500左右
人为入侵在1500~3000左右

pengkun————2016.09.17
***********************************************************************/

#define PI	M_PI
#define NUM_OF_MODEL 2		//定义模板数量
#define nZeroPoint 40		//定义连零点数目
//信号输入
//可控参数设置
/*
// 初步判断可调参数
#define ReduceTail 0.08		// <0.08归为0
#define Judgement_1 1000	// 区域判定1
#define Judgement_2 2000	// 区域判定2
// kmean模版可调参数
#define CosDistance 0.6		// 余弦距离最小值
*/
__para g_Param; 

//信号预处理
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

// 找出连零起始点
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

// 输出判定值
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
//实时数据平方和开根号
double sum_sqrt(double data[], int nn)
{
	double data_r=0.0;
	for (int i = 1; i<nn; i++)		//不确定i的起始
	{
		data_r += pow(data[2 * i + 1], 2.0);
	//	std::cout << data_r<<", "<<data[i]<<". ";
	}
	data_r = sqrt(data_r);
	return data_r;
}
//模板平方和开根号
double model_sqrt(double data[], int nn)
{
	double data_r=0.0;
	for (int i = 1; i<nn; i++)		//不确定i的起始
	{
		data_r += pow(data[i], 2);
	//	std::cout << data_r<<", "<<data[i]<<". ";
	}
	data_r = sqrt(data_r);
	return data_r;
}

//计算余弦距离
double caculate_value(double data[], double data_sqrts,double models[], double model_sqrts, int ModelSize)
{
	double add_result = 0.0;
	double cosdis = 0;
	int i;
	for (i = 1; i<ModelSize; i++)		//不确定i的起始
	{
		add_result += data[2 * i + 1] * models[i];
	}
	cosdis = add_result/(data_sqrts * model_sqrts);    // 128*????????
	//printf("%d,%f,%d\t", add_result, data_sqrts, cosdis);
	return(cosdis);

}
//比较信号与模板之间距离大小
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

//判断信号类型
// adjuestment
int choose_result(int chooses)
{
	int value = 0;
	if ( chooses > 200)	value = 10000;
	else if (chooses > 100)	value = 100;
	else value = 10;
	return value;
}

//信号分类
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
	// 这是用于初步判断
	double* data_pro1;			//用于重构10240个数据，进行事件初步判断
	// 这是用于kmean
	double* data_pro2;			//用于重构10240个数据，使其与算法匹配
	//double DataBuf[N*20];	//输入的10240个数据
	//信号预处?
	//获取数据时域FFT的范围
	int ShowOut;
	NFFT = (int)pow(2.0, ceil(log((double)nDataSize) / log(2.0)));	
	data_pro1 = DataFill(DataBuf, nDataSize, NFFT);	
	DataFFT(data_pro1, NFFT);
	DataReduce(data_pro1, NFFT/2, NFFT/2, g_Param.fReduceTail);//ReduceTail);
	data_end=DataZeroJudge(data_pro1, nZeroPoint, NFFT/2 );//ZeroPoint, NFFT/2 );
	g_Param.DataEnd = data_end;
	delete[] data_pro1;
	preout=GetInterval(data_end, g_Param.nJudgement_1, g_Param.nJudgement_2);	//输出区间代表:100,2,1000

	Param.DataEnd = g_Param.DataEnd;	// 输出展示
	
	if(preout == 2)
	{
		// 进行kmean分析
		data_pro2 = DataFill(DataBuf, nDataSize, NFFT);
		DataFFT(data_pro2, NFFT);
		DataCut(data_pro2, ModelSize, NFFT/2);	//之前代码这里有误，但对数据处理几乎无影响
		//printf("\n\n");
		//计算余弦距离
		out = classify(data_pro2, EnvironTemp, IntrudeTemp, ModelSize, ShowOut);
		g_Param.ShowOutPos1 = ShowOut;
		//判断信号类型
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

