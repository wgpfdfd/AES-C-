#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#define uint unsigned int
using namespace std;

static const unsigned int sixteen_ten[16] = { 0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xa,0xb,0xc,0xd,0xe,0xf };

static const unsigned int graph[16][16] =
{ 0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };   //S box


static const unsigned int inverse_graph[16][16] =
{   0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
	0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
	0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
	0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
	0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
	0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
	0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
	0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
	0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
	0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
	0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
	0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
	0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
	0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
	0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
	0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };  //逆 S box

static const unsigned int rcon[10][8] = {
	0,1,0,0,0,0,0,0,
	0,2,0,0,0,0,0,0,
	0,4,0,0,0,0,0,0,
	0,8,0,0,0,0,0,0,
	1,0,0,0,0,0,0,0,
	2,0,0,0,0,0,0,0,
	4,0,0,0,0,0,0,0,
	8,0,0,0,0,0,0,0,
	1,11,0,0,0,0,0,0,
	3,6,0,0,0,0,0,0 };  //rcon box


	// 域运算辅助函数
unsigned int addDel(unsigned int a, unsigned int b)
{
	unsigned int res = 0x00;
	res = a ^ b;
	return res;
}

unsigned char mod(unsigned char a)
{
	return (a << 1) ^ ((a & 0x80) ? 0x1B : 0);
}
static void TwoOneSymbol(unsigned int* input, unsigned int* result)  //1*31------> 2*16
{
	for (int i = 0; i < 16; i++)
	{
		*result = *input * 16 + *(input + 1);
		if (i != 15)
		{
			result++;
			input += 2;
		}
	}
}

static void OneTwoSymbol(unsigned int* input, unsigned int* result)   //2*16------> 1*31
{
	for (int i = 0; i < 16; i++)
	{
		*result = *input / 16;
		result++;
		*result = *input % 16;
		if (i != 15)
		{
			result++;
			input++;
		}

	}
}
unsigned char multiply(unsigned char a, unsigned char b)
{
	unsigned char res = 0x00;
	unsigned char bitAnd = 0x01;
	unsigned char addList[8] = { a };
	for (int i = 0; i < 7; i++)
	{
		addList[i + 1] = mod(addList[i]);
	}

	for (int i = 0; i < 8; i++)
	{
		if (b & bitAnd)
		{
			res = addDel(res, addList[i]);
		}
		bitAnd = bitAnd << 1;
	}
	return res;
}


static void MatrixToArray(unsigned int* array, unsigned int* result)   //4*4*2 ---->1*31
{
	for (int i = 0; i < 16; i++)
	{
		*result = sixteen_ten[*array / 16];
		result++;
		*result = sixteen_ten[*array % 16];
		if (i != 15)
		{
			result++;
			array++;
		}
	}
}

static void ArrayToMatrix(unsigned int* array, unsigned int result[4][4])   //1*32 ---->  4*4*2
{
	for (int i = 0; i < 16; i++)
	{
		result[i / 4][i % 4] = (*array * 16 + *(array + 1));
		if (i != 15)
		{
			array += 2;
		}

	}
}

static void SubBytes(unsigned int *signal_input,int size)    //字节代换
{
	unsigned int temp[32] = { 0 };
	for (int i = 0; i < size / 2; i++)
	{
		temp[i] = graph[*(signal_input)][*(signal_input + 1)];//  search in graph
		*signal_input = temp[i] / 16;    //取十位
		*(signal_input + 1) = temp[i] % 16;  //取个位
		if (i != size / 2 - 1)
		{
			signal_input += 2;
		}
	}

}

void ShiftRows(unsigned int* signal_input)               //shift rows
{
	unsigned int temp_array[2] = { 0 };
	signal_input += 2;  // start by S1
	for (int i = 0; i < 2; i++)    //second line shift
	{
		temp_array[i] = *signal_input;
		*signal_input = *(signal_input + 8);
		*(signal_input + 8) = *(signal_input + 16);
		*(signal_input + 16) = *(signal_input + 24);
		*(signal_input + 24) = temp_array[i];
		signal_input++;             //when finished ,address of signal_input address = S2
	}

	for (int i = 0; i < 2; i++)    //third line shift
	{
		temp_array[i] = *signal_input;
		*signal_input = *(signal_input + 16);
		*(signal_input + 16) = temp_array[i];
		temp_array[i] = *(signal_input + 8);
		*(signal_input + 8) = *(signal_input + 24);
		*(signal_input + 24) = temp_array[i];
		signal_input++;
	}

	for (int i = 0; i < 2; i++)    //forth line shift
	{
		temp_array[i] = *signal_input;
		*signal_input = *(signal_input + 24);
		*(signal_input + 24) = *(signal_input + 16);
		*(signal_input + 16) = *(signal_input + 8);
		*(signal_input + 8) = temp_array[i];
		if (i != 1)
		{
			signal_input++;
		}
	}

}

void MixColumns(unsigned int* input_signal)   //mix columns
{
	unsigned int tmp[4] = {0};
	unsigned int long_array[32] = {0};
	unsigned int short_array[16] = {0};
	int i, j;
	for (int k = 0; k < 32; k++)
	{
		long_array[k] = *input_signal;
		if (k != 31)
		{
			input_signal++;
		}
	}
	input_signal -= 31;
	TwoOneSymbol(long_array, short_array);        // 1*32--->2*16
	for (j = 0; j < 4; j++)      //按列处理
	{

		for (i = 0; i < 4; i++) {

			tmp[i] = short_array[j * 4 + i];      //每一列中的每一个字节拷贝到tmp中
		}
		for (i = 0; i < 4; i++) {

			short_array[j * 4 + i] = multiply(0x02, tmp[i])      //矩阵计算，加法为异或
				^ multiply(0x03, tmp[(i + 1) % 4])
				^ multiply(0x01, tmp[(i + 2) % 4])
				^ multiply(0x01, tmp[(i + 3) % 4]);
		}
	}
	OneTwoSymbol(short_array, long_array);          //2*16--->1*32
	for (int m = 0; m < 32; m++)
	{
		*input_signal = long_array[m];
		if (m != 31)
		{
			input_signal++;
		}
	}
}


static void str_left_shift(unsigned int* array)  //左移2位  used in key expand
{
	unsigned int temp[2] = { 0 };
	temp[0] = *array;
	temp[1] = *(array + 1);
	for (int i = 0; i < 6; i++)
	{
		*array = *(array + 2);
		array++;
	}
	*array = temp[0];
	*(array + 1) = temp[1];
}


void key_array_expand(unsigned int* initial_key_array, unsigned int* key_addition)     //密钥拓展
{
	unsigned int temp_key_array[8] = { 0 };
	unsigned int temp_key_addition[352] = { 0 };
	for (int i = 0; i < 32; i++)    //first 4( W[0]--W[4])  for W[i] = K[i],K[i+1],K[i+2],K[i+3]     need  10------>16
	{
		temp_key_addition[i] = *initial_key_array;
		initial_key_array++;
	}

	for (int j = 4; j < 44; j++)
	{
		if (j % 4 == 0)
		{
			for (int m = 0; m < 8; m++)
			{
				temp_key_array[m] = temp_key_addition[8 * (j - 1) + m];   //temp_key_array = W[i-1]
			}

			str_left_shift(temp_key_array); //左移两位
			SubBytes(temp_key_array,8);  //replace in S box

			for (int s = 0; s < 8; s++)
			{
				temp_key_array[s] = temp_key_array[s] ^ rcon[j / 4 - 1][s];   //xor with rcon matrix
			}

			for (int n = 0; n < 8; n++)
			{
				temp_key_addition[8 * j + n] = temp_key_addition[8 * (j - 4) + n] ^ temp_key_array[n];   // W[i-4] xor T( W[i-1] )
			}


		}

		else    //W[j] = W[j-4] xor W[j-1] 
		{
			for (int k = 0; k < 8; k++)
			{
				temp_key_addition[8 * j + k] = temp_key_addition[8 * (j - 4) + k] ^ temp_key_addition[8 * (j - 1) + k];   //W[i-4] xor W[i-1]
			}
		}

	}

	for (int p = 0; p < 352; p++)  //put result into key_addition
	{
		*key_addition = temp_key_addition[p];
		key_addition++;
	}
}


static void InvSubBytes(unsigned int* signal_input, int size)    //逆字节代换
{
	unsigned int temp[32] = { 0 };
	for (int i = 0; i < size / 2; i++)
	{
		temp[i] = inverse_graph[*(signal_input)][*(signal_input + 1)];   //search in graph
		*signal_input = temp[i] / 16;
		*(signal_input + 1) = temp[i] % 16;
		if (i != size / 2 - 1)
		{
			signal_input += 2;
		}
	}

}

void InvShiftRows(unsigned int* signal_input)  //逆行位移
{
	unsigned int temp_array[2];
	signal_input += 2;   // start by S1
	for (int i = 0; i < 2; i++)    //second line shift
	{
		temp_array[i] = *signal_input;
		*signal_input = *(signal_input + 24);
		*(signal_input + 24) = *(signal_input + 16);
		*(signal_input + 16) = *(signal_input + 8);
		*(signal_input + 8) = temp_array[i];
		signal_input++;             //when finished,  address of signal_input = S2
	}

	for (int i = 0; i < 2; i++)    //third line shift
	{
		temp_array[i] = *signal_input;
		*signal_input = *(signal_input + 16);
		*(signal_input + 16) = temp_array[i];
		temp_array[i] = *(signal_input + 8);
		*(signal_input + 8) = *(signal_input + 24);
		*(signal_input + 24) = temp_array[i];
		signal_input++;
	}

	for (int i = 0; i < 2; i++)    //forth line shift

	{
		temp_array[i] = *signal_input;
		*signal_input = *(signal_input + 8);
		*(signal_input + 8) = *(signal_input + 16);
		*(signal_input + 16) = *(signal_input + 24);
		*(signal_input + 24) = temp_array[i];
		if (i != 1)
		{
			signal_input++;
		}
	}
}

void InvMixColumns(unsigned int* input_signal)  //逆列混合
{
	unsigned int tmp[4] = {0};
	int i, j;
	unsigned int long_array[32] = { 0 };
	unsigned int short_array[16] = { 0 };
	for (int k = 0; k < 32; k++)
	{
		long_array[k] = *input_signal;
		if (k != 31)
		{
			input_signal++;
		}
	}
	input_signal -= 31;
	TwoOneSymbol(long_array, short_array);    // 1*32--->2*16
	for (j = 0; j < 4; j++)      //按列处理
	{

		for (i = 0; i < 4; i++) {

			tmp[i] = short_array[j * 4 + i];      //每一列中的每一个字节拷贝到tmp中
		}
		for (i = 0; i < 4; i++) {

			short_array[j * 4 + i] = multiply(0x0e, tmp[i])      //矩阵计算，加法为异或
				^ multiply(0x0b, tmp[(i + 1) % 4])
				^ multiply(0x0d, tmp[(i + 2) % 4])
				^ multiply(0x09, tmp[(i + 3) % 4]);
		}
	}
	OneTwoSymbol(short_array, long_array);   //2*16--->1*32
	for (int m = 0; m < 32; m++)
	{
		*input_signal = long_array[m];
		if (m != 31)
		{
			input_signal++;
		}
	}
}

void AddRoundKey(unsigned int* input_signal, unsigned int* key_addition)   //轮密钥加
{
	for (int i = 0; i < 32; i++)
	{
		*input_signal = *input_signal ^ *key_addition;   //xor
		if (i != 31)
		{
			input_signal++;
			key_addition++;
		}

	}
}

static void PrintfKeyAddition(unsigned int* keyaddition)   //key printf function
{
	unsigned int temp[352];
	unsigned int printfkeyaddition[4][4][11] = {0};

	for (int m = 0; m < 352; m++)      //put key in temp
	{
		temp[m] = *keyaddition;
		if (m != 351)
		{
			keyaddition++;
		}
	}
	for (int k = 0; k < 11; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				printfkeyaddition[i][j][k] = temp[32 * k + 2 * i + j * 8] * 16 + temp[32 * k + 2 * i + j * 8 + 1];
				printf("%x,", printfkeyaddition[i][j][k]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
}


static void printfinmatrix(unsigned int* signal)    //printf function  1*32----4*4*2
{
	unsigned int temp[32] = {0};
	unsigned int printfsignal[4][4] = { 0 };
	for (int m = 0; m < 32; m++)
	{
		temp[m] = *signal;
		if (m != 31)
		{
			signal++;
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			printfsignal[i] [j] = temp[2 * i + j * 8] * 16 + temp[2 * i + j * 8 + 1];
			printf("%x,", printfsignal[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}


int main()
{
	unsigned int input_signal[16] = {
		0x32,0x43,0xf6,0xa8,
		0x88,0x5a,0x30,0x8d,
		0x31,0x31,0x98,0xa2,
		0xe0,0x37,0x07,0x34 };             //S0=87, S1=F2,S2=4D.....

	unsigned int key_array[16] = {
		0x2b,0x7e,0x15,0x16,
		0x28,0xae,0xd2,0xa6,
		0xab,0xf7,0x15,0x88,
		0x09,0xcf,0x4f,0x3c };
	unsigned int temp_output[16];
	unsigned int long_array[32] = { 0 };     //key  W[0]--W[3]
	unsigned int key_addition[352] = { 0 };//key after expand
	unsigned int long_signal[32] = { 0 };  //signal in one array
	MatrixToArray(key_array, long_array);
	MatrixToArray(input_signal, long_signal);


	printf("\ninput_signal:\n");
	printfinmatrix(long_signal);
	key_array_expand(long_array, key_addition);  //密钥拓展

	unsigned int temp_key[32] = { 0 };
	printf("\nkey after expand:\n");
	PrintfKeyAddition(key_addition);


	//encropy processing
	printf("\nencropy starting...\n\n");
	AddRoundKey(long_signal, long_array);
	printfinmatrix(long_array);
	printf("\nround 1 start by matrix:\n");
	printfinmatrix(long_signal);
	for (int i = 1; i < 10; i++)         //round  1--9
	{
		printf("\nin %d round:\n\n", i);
		SubBytes(long_signal,32);
		printf("\nafter subBytes:\n");
		printfinmatrix(long_signal);
		printf("\nafter ShiftRows:\n");
		ShiftRows(long_signal);
		printfinmatrix(long_signal);
		printf("\nafter MixColumns:\n");
		MixColumns(long_signal);
		printfinmatrix(long_signal);
		for (int j = 0; j < 32; j++)
		{
			temp_key[j] = key_addition[i * 32 + j];     //for addRoundKey    W[i]---W[i+3]
		}
		AddRoundKey(long_signal, temp_key);
		printfinmatrix(long_signal);
	}
	for (int i = 0; i < 32; i++)
	{
		temp_key[i] = key_addition[320 + i];     //W[40]--W[43]
	}
	SubBytes(long_signal, 32);
	ShiftRows(long_signal);
	AddRoundKey(long_signal, temp_key);        //last round
	printf("\nencropy array = :\n");
	printfinmatrix(long_signal);


	// dencropy processing
	printf("\ndencropy starting...\n\n");
	printf("\nround 1 start by matrix:\n");
	AddRoundKey(long_signal, temp_key);
	printfinmatrix(long_signal);

	for (int i = 9; i > 0; i--)                 //round  1-- 9
	{
		printf("in %d round:", 10-i);
		printf("\nafter InvShiftRows:\n");
		InvShiftRows(long_signal);
		printfinmatrix(long_signal);
		printf("\nafter SubBytes:\n");
		InvSubBytes(long_signal, 32);
		printfinmatrix(long_signal);
		for (int j = 0; j < 32; j++)
		{
			temp_key[j] = key_addition[i * 32 + j];     //for addRoundKey    W[i]---W[i+3]
		}
		AddRoundKey(long_signal, temp_key);
		printf("\nafter AddRoundKey:\n");
		printfinmatrix(long_signal);
		InvMixColumns(long_signal);
		printf("\nafter InvMixColumns:\n");
		printfinmatrix(long_signal);
	}


	InvShiftRows(long_signal);
	InvSubBytes(long_signal, 32);
	AddRoundKey(long_signal, long_array);      ////last round

	printf("\ndencropy array = :\n");
	printfinmatrix(long_signal);

	return 0;
}


