/*
 * PaillierCrypto.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: PJS
 */


#include "PaillierCrypto.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <wchar.h>
#include <cstring>
#include "SocketException.h"


PaillierCrypto::PaillierCrypto() : mPubKey(NULL)
{
	try {
		ServerSocket Server(3000);
		mSSocket = new ServerSocket();
		Server.accept(*mSSocket);
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	gmp_randinit_default(state);
	mPrvKey = NULL;									// for debugging
}

PaillierCrypto::~PaillierCrypto()
{
	paillier_freepubkey(mPubKey);
	paillier_freeprvkey(mPrvKey);					// for debugging
	delete(mSSocket);
}

paillier_pubkey_t*			PaillierCrypto::GetPubKey()		{return mPubKey;}

void PaillierCrypto::PreComputation(pre_t *tPre, skle_t* tSkle, out_t *tOut, in_t* tIn, unsigned short idx)
{
	paillier_plaintext_t pT;		// 임시변수

	#ifdef _DEBUG_INIT_1
	mpz_inits(pT.m, tPre->pN1.m, tPre->pN2.m, tPre->pL.m, tPre->c1.c, tPre->cN1.c, tPre->cNk.c, NULL);
	#else
	mpz_init2(pT.m, GMP_N_SIZE+1);
	mpz_init2(tPre->pN1.m, GMP_N_SIZE+1);
	mpz_init2(tPre->pN2.m, GMP_N_SIZE+1);
	mpz_init2(tPre->pL.m, GMP_N_SIZE+1);
	mpz_init2(tPre->c1.c, 3*GMP_N_SIZE);
	mpz_init2(tPre->cN1.c, 3*GMP_N_SIZE);
	mpz_init2(tPre->cNk.c, 3*GMP_N_SIZE);
	#endif

	for (int i=0 ; i<CLASS_SIZE ; i++)
		#ifdef _DEBUG_INIT_1
		mpz_init(tPre->cCls[i].c);
		#else
		mpz_init2(tPre->cCls[i].c, 2*GMP_N_SIZE*2);
		#endif

	// initialize cS = E(0)
	for (int i=0 ; i<DATA_SIZE ; i++) {
		mpz_set_ui(tOut->cS[i].c, 1);
	}

	// N-1
	mpz_sub_ui(tPre->pN1.m, mPubKey->n, 1);			// mod 연산 추가
	// N-2
	mpz_sub_ui(tPre->pN2.m, mPubKey->n, 2);			// mod 연산 추가

	#ifdef _DEBUG_Initialization
	DebugOut("N\t\t\t", mPubKey->n, idx);
	DebugOut("N-1\t\t\t", tPre->pN1.m, idx);
	DebugOut("N-2\t\t\t", tPre->pN2.m, idx);
	#endif

	// l = 2^(-1) mod N
	mpz_set_ui(pT.m, 2);			// 임시 변수로서 pN1 이용함.
	assert(mpz_invert(tPre->pL.m, pT.m, mPubKey->n));

	#ifdef _DEBUG_Initialization
	DebugOut("l=2^(-1) mod N\t\t", tPre->pL.m, idx);
	#endif

	// c1 = E(1)
	mpz_set_ui(pT.m, 1);													// pT=1
	paillier_enc(&(tPre->c1), mPubKey, &pT, paillier_get_rand_devurandom);

	#ifdef _DEBUG_Initialization
	//DebugOut("c1 = E(1)\t\t", tPre->c1, idx);
	DebugDec("c1 = E(1)\t\t", &(tPre->c1), idx);
	#endif

	// compute (-c_j)' = E(N-c_j) for j=1~v
	for (int j=0 ; j<CLASS_SIZE ; j++) {
		mpz_sub_ui(pT.m, mPubKey->n, tIn->Class[j]);
		//tCfd->cClass_[CLASS_SIZE];
		paillier_enc(&(tPre->cCls[j]), mPubKey, &pT, paillier_get_rand_devurandom);

		#ifdef _DEBUG_Initialization
		//DebugOut("(-c_j)' = E(N-c_j)\t", tPre->cCls[j], idx);
		DebugDec("(-c_j)' = E(N-c_j)\t", &(tPre->cCls[j]), idx);
		#endif
	}

	// cN1 = E(N-1)
	paillier_enc(&(tPre->cN1), mPubKey, &(tPre->pN1), paillier_get_rand_devurandom);

	// cNk = E(N-k)
	mpz_sub_ui(pT.m, mPubKey->n, PARAM_K);							// pT=N-k
	paillier_enc(&(tPre->cNk), mPubKey, &pT, paillier_get_rand_devurandom);

	// Comparison 함수의 input parameter
	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
		#ifdef _DEBUG_INIT_1
		mpz_inits(tSkle->cCntbit[0][i].c, tSkle->cCntbit[1][i].c, NULL);
		#else
		{ mpz_init2(tSkle->cCntbit[0][i].c, 2*GMP_N_SIZE*2);
		  mpz_init2(tSkle->cCntbit[1][i].c, 2*GMP_N_SIZE*2);	}
		#endif
	}

	for (int i=0 ; i<NUM_THREAD_DATA ; i++)
		#ifdef _DEBUG_INIT_1
		//mpz_inits(tSkle->cTRes[i].c, tSkle->cIRes[i].c, tSkle->cCan[i].c, NULL);
		mpz_inits(tSkle->cC[i].c, tSkle->cK[i].c, NULL);
		#else
		{ mpz_init2(tSkle->cTRes[i].c, 2*GMP_N_SIZE*2);
		  mpz_init2(tSkle->cIRes[i].c, 2*GMP_N_SIZE*2);
		  mpz_init2(tSkle->cCan[i].c, 2*GMP_N_SIZE*2);  }
		#endif

	// Cmp func에서 SBD func을 이용하는지 여부를 위한 변수 (data의 SBD -> Cmp의 SBD -> frequency의 SBD -> Cmp의 SBD )
	tSkle->bSkle = false;				// SBD

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(		  pT.m->_mp_alloc ==  GMP_N_SIZE+1);
	assert(tPre->pN1.m->_mp_alloc ==  GMP_N_SIZE+1);
	assert(tPre->pN2.m->_mp_alloc ==  GMP_N_SIZE+1);
	assert(tPre-> pL.m->_mp_alloc ==  GMP_N_SIZE+1);
	assert(tPre-> c1.c->_mp_alloc ==3*GMP_N_SIZE);
	assert(tPre->cN1.c->_mp_alloc ==2*GMP_N_SIZE*2);
	assert(tPre->cNk.c->_mp_alloc ==2*GMP_N_SIZE*2);
	#endif

	return ;
}

void PaillierCrypto::SetPubKey()
{
	char cRecvBuf[KEY_SIZE*2+1]={0,};
	int iRecvSize;

	// distribute public key
	try {
		iRecvSize = mSSocket->recv(cRecvBuf, KEY_SIZE*2);		// string으로 전송
    }
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	#ifdef _DEBUG_Assert
	assert(iRecvSize == KEY_SIZE*2);
	#endif

	#ifdef _DEBUG_Initialization
    std::cout << "[DH] Public Key (Hex)\t\t : " << cRecvBuf << std::endl;
	#endif

	mPubKey = paillier_pubkey_from_hex(cRecvBuf);
	#ifdef _DEBUG_Initialization
    DebugOut("Public Key (Copied)\t", mPubKey->n, 999);
	#endif

	return;
}

// for debugging
void PaillierCrypto::SetPrvKey()
{
	char cRecvBuf[KEY_SIZE*2+1]={0,};
	int iRecvSize;

	// distribute public key
	try {
		iRecvSize = mSSocket->recv(cRecvBuf, KEY_SIZE*2);			// string으로 전송
    }
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	#ifdef _DEBUG_Assert
	assert(iRecvSize == KEY_SIZE*2);
	#endif

	#ifdef _DEBUG_Initialization
    std::cout << "[DH] Private Key (Hex)\t\t : " << cRecvBuf << std::endl;
	#endif

	mPrvKey = paillier_prvkey_from_hex(cRecvBuf, mPubKey);
	#ifdef _DEBUG_Initialization
    DebugOut("Private Key (Copied)\t", mPrvKey->lambda, 999);
	#endif

	return ;
}

bool PaillierCrypto::inputDatasetQuery(in_t *tIn)
{
	std::fstream DatasetFile("dataset.dat", std::fstream::in);
	std::fstream InqueryFile("inquery.dat", std::fstream::in);

	#ifdef _DEBUG_Assert
	assert(DatasetFile.is_open());
	assert(InqueryFile.is_open());
	#endif

	std::string stTmp;
	paillier_plaintext_t pTmp;
	std::string stToken;

	#ifdef _DEBUG_INIT_1
	mpz_init(pTmp.m);
	#else
	mpz_init2(pTmp.m, 1);
	#endif

	printf("\n************   Initialization   ************\n\n");

	// input dataset
	for (int i=0 ; i<DATA_SIZE ; i++) {
		std::getline(DatasetFile, stTmp);
		#ifdef _DEBUG_Initialization
		std::cout << "[DH] Dataset:\t\t" << stTmp << std::endl;
		#endif

		std::stringstream stStream(stTmp);
		for (int j=0 ; j<DATA_DIMN ; j++) {
			stStream >> stToken ;
			mpz_set_ui(pTmp.m, std::atoi(stToken.c_str()));
			paillier_enc(&(tIn->cData[i][j]), mPubKey, &pTmp, paillier_get_rand_devurandom);

			#ifdef _DEBUG_Initialization
			std::cout << "[DH] Dataset:\t\t" << stToken << std::endl;
			gmp_printf("[DH] Dataset:\t\t%ZX\n", &pTmp);
			gmp_printf("[DH] Encrypted Dataset: %ZX\n", &(tIn->cData[i][j]));
			#endif
		}
		stStream >> stToken ;
		mpz_set_ui(pTmp.m, std::atoi(stToken.c_str()));
		paillier_enc(&(tIn->cClass[i]), mPubKey, &pTmp, paillier_get_rand_devurandom);

		#ifdef _DEBUG_Initialization
		std::cout << "[DH] Dataset Class:\t" << stToken << std::endl;
		gmp_printf("[DH] Dataset Class:\t%ZX\n", &pTmp);
		gmp_printf("[DH] Encrypted Dataset Class: %ZX\n", &(tIn->cClass[i]));
		#endif

		#ifdef _DEBUG_Assert
		assert(i<=DATA_SIZE);
		#endif
	}

	// input query
	std::getline(InqueryFile, stTmp);
	#ifdef _DEBUG_Initialization
	std::cout << "[DH] Query:\t\t" << stTmp << std::endl;
	#endif

	std::stringstream stStream(stTmp);
	for (int j=0 ; j<DATA_DIMN ; j++) {
		stStream >> stToken ;
		mpz_set_ui(pTmp.m, std::atoi(stToken.c_str()));
		paillier_enc(&(tIn->cQuery[j]), mPubKey, &pTmp, paillier_get_rand_devurandom);

		#ifdef _DEBUG_Initialization
		std::cout << "[DH] Query:\t\t" << stToken << std::endl;
		gmp_printf("[DH] Query:\t\t%ZX\n", &pTmp);
		gmp_printf("[DH] Encrypted Query:   %ZX\n", &(tIn->cQuery[j]));
		#endif
	}

	DatasetFile.close();
	InqueryFile.close();

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(pTmp.m->_mp_alloc == 1);
	#endif

	return true;
}

//void PaillierCrypto::SquaredDist(out_t* tOut, in_t *tIn, pre_t* tPre, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::SquaredDist(out_t* tOut, in_t *tIn, pre_t* tPre, unsigned short idx, thp_t* tRecv)
{
	paillier_ciphertext_t cTmp1, cTmp2, cQN1[DATA_DIMN];

	#ifdef _DEBUG_INIT_1
	mpz_inits(cTmp1.c, cTmp2.c, NULL);
	#else
	mpz_init2(cTmp1.c, 2*GMP_N_SIZE*2);
	mpz_init2(cTmp2.c, 2*GMP_N_SIZE*2);
	#endif

	for (int j=0 ; j<DATA_DIMN ; j++)
		#ifdef _DEBUG_INIT_1
		mpz_init(cQN1[j].c);
		#else
		mpz_init2(cQN1[j].c, 2*GMP_N_SIZE);
		#endif

	// N-1
	#ifdef _DEBUG_SquaredDist
	DebugOut("N   \t", mPubKey->n, idx);	DebugOut("N-1 \t", tPre->pN1.m, idx);
	#endif

	// q'_j ^(N-1)  for j=1~m
	for (int j=0 ; j<DATA_DIMN ; j++) {
		paillier_exp(mPubKey, &cQN1[j], &(tIn->cQuery[j]), &(tPre->pN1));
		#ifdef _DEBUG_SquaredDist
		//DebugOut("E(-q[j])", cQN1[j].c, idx);
		DebugDec("-q[j]\t", &cQN1[j], idx);
		#endif
	}

	for (int i=tOut->iRan[idx] ; i<tOut->iRan[idx+1] ; i++) {
		// initialize s' = E(0)
		#ifdef _DEBUG_SquaredDist
		//DebugOut("Encrypted s=0 (initialization)\t", tOut->cS[i].c, idx);
		DebugDec("s=0 (initialization)\t\t", &(tOut->cS[i]), idx);
		#endif

		for (int j=0 ; j<DATA_DIMN ; j++) {
			// tmp1' = d'_(i,j)*q'_j^(N-1) = E(d_(i,j)-q_j) for j=1~m, i=1~n
			paillier_mul(mPubKey, &cTmp1, &(tIn->cData[i][j]), &cQN1[j]);
			#ifdef _DEBUG_SquaredDist
			//DebugOut("Encrypted tmp1 = E(d_i,j-q_j)", cTmp1.c, idx);
			DebugDec("tmp1 = d_i,j-q_j\t\t\t", &cTmp1, idx);
			#endif

			// tmp2' = (tmp1'_j)^2 =  E((d_(i,j)-q_j)^2) for j=1~m, i=1~n
			//SecMul(&cTmp2, &cTmp1, idx, tRecv, tSend, ucSendPtr);
			SecMul(&cTmp2, &cTmp1, idx, tRecv);
			#ifdef _DEBUG_SquaredDist
			//DebugOut("Encrypted tmp2 = E((d-q)^2)", cTmp2.c, idx);
			DebugDec("tmp2 = (d-q)^2\t\t\t\t", &cTmp2, idx);
			#endif

			// tmp3' =  E((d_(i,1)-q_1)^2) + E((d_(i,2)-q_2)^2) + ... + E((d_(i,m)-q_m)^2)
			paillier_mul(mPubKey, &(tOut->cS[i]), &(tOut->cS[i]), &cTmp2);
			#ifdef _DEBUG_SquaredDist
			//DebugOut("Encrypted s = E((d-q)^2 + ...)", tOut->cS[i].c, idx);
			DebugDec("s = (d-q)^2 + ...\t\t\t", &(tOut->cS[i]), idx);
			#endif
		}

		#ifdef _DEBUG_SquaredDist
		//DebugOut("Encrypted Squared Distance s = E((d-q)^2 + ...)", tOut->cS[i].c, idx);
		DebugDec("Encrypted Squared Distance s = (d-q)^2 + ...\t\t\t", &(tOut->cS[i]), idx);
		#endif

	}

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(  cTmp1.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(  cTmp2.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(cQN1[0].c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

	return;
}

//void PaillierCrypto::SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, paillier_ciphertext_t* cEb, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, paillier_ciphertext_t* cEb, unsigned short idx, thp_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+2*ENC_SIZE]={0,};
	unsigned char bA[ENC_SIZE]={0,}, bB[ENC_SIZE]={0,};
	short sLen;

	paillier_plaintext_t  pRa, pRb, pNt;
	paillier_ciphertext_t cA, cB, cT, cH;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pRa.m, pRb.m, pNt.m, cA.c, cB.c, cT.c, cH.c, NULL);
	#else
	mpz_init2(pRa.m, GMP_N_SIZE);
	mpz_init2(pRb.m, GMP_N_SIZE);
	mpz_init2(pNt.m, 2*GMP_N_SIZE);
	mpz_init2(cA.c, 2*GMP_N_SIZE*2);
	mpz_init2(cB.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cH.c, 2*GMP_N_SIZE);
	#endif

	// parameter : E(a), E(b)
	#ifdef _DEBUG_SecureMultiplication
	DebugDec("input a\t\t\t", cEa, idx);
	DebugDec("input b\t\t\t", cEb, idx);
	#endif

	// 0 <= r_a, r_b < N  (Pick	two random numbers)
	mpz_urandomm(pRa.m, state, mPubKey->n);
	mpz_urandomm(pRb.m, state, mPubKey->n);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("r_a (random number 1)\t\t", pRa.m, idx);
	DebugOut("r_b (random number 2)\t\t", pRb.m, idx);
	#endif

	// E(r_a), E(r_b)
	paillier_enc(&cA, mPubKey, &pRa, paillier_get_rand_devurandom);
	paillier_enc(&cB, mPubKey, &pRb, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("E(r_a)\t\t\t\t", cA.c, idx);
	DebugOut("E(r_b)\t\t\t\t", cB.c, idx);
	#endif

	// a' = E(a)*E(r_a) = E(a+r_a)
	// b' = E(b)*E(r_b) = E(b+r_b)
	paillier_mul(mPubKey, &cA, &cA, cEa);
	paillier_mul(mPubKey, &cB, &cB, cEb);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("a'=E(a+r_a)\t\t\t", cA.c, idx);
	DebugOut("b'=E(b+r_b)\t\t\t", cB.c, idx);
	DebugDec("a + r_a\t\t\t", &cA, idx);
	DebugDec("b + r_b\t\t\t", &cB, idx);
	#endif

	// send a', b'
	paillier_ciphertext_to_bytes(bA, ENC_SIZE, &cA);
	paillier_ciphertext_to_bytes(bB, ENC_SIZE, &cB);
	SetSendMsg(ucSendPtr, bA, bB, idx, COM_MUL2, ENC_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SecureMultiplication
    DebugCom("(buf) a'=E(a+r_a), b'=E(b+r_b)\t", ucSendPtr, 2*ENC_SIZE+HED_SIZE, idx);
	#endif

//    mSendMtx->lock();
//    mSendQueue->push(ucSendPtr);
//    mSendMtx->unlock();
//    mSendCV->notify_all();

	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// t1 = E(a)^(N-r_b) = -r_b * E(a)
	// ==> Nt = N-r_b, A = E(a)^Nt
	mpz_sub(pNt.m, mPubKey->n, pRb.m);
	paillier_exp(mPubKey, &cA, cEa, &pNt);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("t1 = E(a)^(N-r_b) = -r_b * E(a)", cA.c, idx);
	DebugDec("t1 = E(a)^(N-r_b) = -r_b * E(a)", &cA, idx);
	#endif

	// t2 = E(b)^(N-r_a) = -r_a * E(b)
	// ==> Nt = N-r_a, B = E(b)^Nt
	mpz_sub(pNt.m, mPubKey->n, pRa.m);
	paillier_exp(mPubKey, &cB, cEb, &pNt);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("t2 = E(b)^(N-r_a) = -r_a * E(b)", cB.c, idx);
	DebugDec("t2 = E(b)^(N-r_a) = -r_a * E(b)", &cB, idx);
	#endif

	// tmp_Res1 = E(-r_b*a) * E(-r_a*b) = E(- r_b*a - r_a*b)
	// ==> tmp_Res1 = t1 * t2
    paillier_mul(mPubKey, cRes, &cA, &cB);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("tmp_Res1 = t1 * t2\t\t", cRes->c, idx);
	DebugDec("tmp_Res1 = t1 * t2\t\t", cRes, idx);
	#endif

	// E(r_a*r_b)^(N-1) = E(- r_a*r_b) = E( (N-r_a)*r_b )
	// ==> Nt = N-r_a, Nt = Nt * r_b,
	// ==> T = E(Nt)
	mpz_sub(pNt.m, mPubKey->n, pRa.m);
	mpz_mul(pNt.m, pNt.m, pRb.m);
	mpz_mod(pNt.m, pNt.m, mPubKey->n);
	paillier_enc(&cT, mPubKey, &pNt, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("Tmp = E(r_a*r_b)^(N-1)\t\t", cT.c, idx);
	DebugDec("Tmp = E(r_a*r_b)^(N-1)\t\t", &cT, idx);
	#endif

	// tmp_Res2 =tmp_Res1 * T
	paillier_mul(mPubKey, cRes, cRes, &cT);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("E(-r_b*a)*E(-r_a*b)*E(-r_a*r_b)", cRes->c, idx);
	DebugDec("E(-r_b*a)*E(-r_a*b)*E(-r_a*r_b)", cRes, idx);
	#endif

	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// E( (a+r_a)*(b+r_b) ) = E( ab + a*r_b + b*r_a + r_a*r_b )
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SecureMultiplication
    DebugCom("Encrypted h' (Hex)\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	// a' = E(a+r_a)
	#ifdef _DEBUG_Assert
    assert((*(ucRecvPtr+2)==COM_MUL2)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(&cH, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("Encrypted h' (Copied)\t\t", cH.c, idx);
	#endif

	// E(ab) = E( ab + a*r_b + b*r_a + r_a*r_b ) * E( - a*r_b - b*r_a - r_a*r_b )
    paillier_mul(mPubKey, cRes, cRes, &cH);
	#ifdef _DEBUG_SecureMultiplication
	//DebugOut("E(a*b)", cRes->c, idx);
	DebugDec("E(a*b)\t\t\t\t", cRes, idx);
	#endif

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(pRa.m->_mp_alloc ==   GMP_N_SIZE);
	assert(pRb.m->_mp_alloc ==   GMP_N_SIZE);
	assert(pNt.m->_mp_alloc == 2*GMP_N_SIZE);
	assert( cA.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert( cB.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert( cT.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert( cH.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;
}

//void PaillierCrypto::SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::SecMul(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cEa, unsigned short idx, thp_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bA[ENC_SIZE]={0,};
	short sLen;

	paillier_plaintext_t  pRa, pNt;
	paillier_ciphertext_t cA, cT, cH;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pRa.m, pNt.m, cA.c, cT.c, cH.c, NULL);
	#else
	mpz_init2(pRa.m, GMP_N_SIZE);
	mpz_init2(pNt.m, GMP_N_SIZE+1);
	mpz_init2(cA.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cH.c, 2*GMP_N_SIZE);
	#endif

	// parameter : E(a)
	#ifdef _DEBUG_SecureMultiplication
	DebugDec("input a\t\t\t", cEa, idx);
	#endif

	// 0 <= r_a < N  (Pick	a random number)
	mpz_urandomm(pRa.m, state, mPubKey->n);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("r (random number)\t\t\t", pRa.m, idx);
	#endif

	// E(r_a)
	paillier_enc(&cA, mPubKey, &pRa, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("E(r_a)\t\t\t\t\t", cA.c, idx);
	#endif

	// a' = E(a)*E(r_a) = E(a+r_a)
	paillier_mul(mPubKey, &cA, &cA, cEa);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("a'=E(a+r_a)\t\t\t\t", cA.c, idx);
	DebugDec("a + r_a\t\t\t\t", &cA, idx);
	#endif

	// send a'
	paillier_ciphertext_to_bytes(bA, ENC_SIZE, &cA);
	SetSendMsg(ucSendPtr, bA, idx, COM_MUL1, ENC_SIZE);
	#ifdef _DEBUG_SecureMultiplication
    DebugCom("(buf) a'=E(a+r_a) (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

//    mSendMtx->lock();
//    mSendQueue->push(ucSendPtr);
//    mSendMtx->unlock();
//    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?

	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// t1 = E(a)^(N-2r_a) = -2r_a * E(a)
	// ==> Nt = 2r_a (mod N), Nt = N - 2r_a, A ;= E(a)^Nt
	mpz_mul_2exp(pNt.m, pRa.m, 1);
	mpz_mod(pNt.m, pNt.m, mPubKey->n);
	mpz_sub(pNt.m, mPubKey->n, pNt.m);
	paillier_exp(mPubKey, &cA, cEa, &pNt);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("t1 = E(a)^(N-2r_a) = -2r_a * E(a)\t", cA.c, idx);
	DebugDec("t1 = E(a)^(N-2r_a) = -2r_a * E(a)\t", &cA, idx);
	#endif

	// t2 = E(-r_a^2) = E(N - r_a^2)
	// ==> Nt = r_a^2 (mod N), Nt = N - r_a^2
	// ==> T = E(Nt)
	mpz_powm_ui(pNt.m, pRa.m, 2, mPubKey->n);
	mpz_sub(pNt.m, mPubKey->n, pNt.m);
	paillier_enc(&cT, mPubKey, &pNt, paillier_get_rand_devurandom);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("t2 = E(-r_a^2)\t\t\t\t", cT.c, idx);
	DebugDec("t2 = E(-r_a^2)\t\t\t\t", &cT, idx);
	#endif

	// tmp_Res2 = T * A
	paillier_mul(mPubKey, cRes, &cT, &cA);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("E(-2r_a) * E(-r_a^2)\t\t\t", cRes->c, idx);
	DebugDec("E(-2r_a) * E(-r_a^2)\t\t\t", cRes, idx);
	#endif

	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// E( (a+r_a)^2 ) = E( a^2 + 2*a*r_a + r_a^2 )
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SecureMultiplication
    DebugCom("Encrypted h' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	// a' = E(a+r_a)
	#ifdef _DEBUG_Assert
    assert((*(ucRecvPtr+2)==COM_MUL1)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(&cH, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_SecureMultiplication
	DebugOut("Encrypted h' (Copied)\t\t\t", cH.c, idx);
	#endif

	// E(ab) = E( a^2 + 2*a*r_a + r_a^2 ) * E(-2r_a*a - r_a^2)
    paillier_mul(mPubKey, cRes, cRes, &cH);
	#ifdef _DEBUG_SecureMultiplication
	//DebugOut("E(a^2)", cRes->c, idx);
	DebugDec("E(a^2)\t\t\t\t\t", cRes, idx);
	#endif

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
	assert(pRa.m->_mp_alloc == GMP_N_SIZE);
	assert(pNt.m->_mp_alloc == GMP_N_SIZE+1);
	assert( cA.c->_mp_alloc ==2*GMP_N_SIZE*2);
	assert( cT.c->_mp_alloc ==2*GMP_N_SIZE*2);
	assert(	cH.c->_mp_alloc ==2*GMP_N_SIZE);
	#endif

	memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;
}

//void PaillierCrypto::SecBitDecomp(out_t* tOut, skle_t* tSkle, pre_t* tPre, bool bFirst, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::SecBitDecomp(out_t* tOut, skle_t* tSkle, pre_t* tPre, bool bFirst, unsigned short idx, thp_t* tRecv)
{
	paillier_ciphertext_t cT, cZ;
	paillier_ciphertext_t *cEX;
	paillier_ciphertext_t *cEXi;
	int iGamma;
	int iSrt, iEnd;
	int iLen;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cT.c, cZ.c, NULL);
	#else
	mpz_init2(cT.c, 2*GMP_N_SIZE);
	mpz_init2(cZ.c, 2*GMP_N_SIZE*2);
	#endif

	if (tSkle->bSkle) {
		cEX		= &(tOut->cCnt);	//tSkle->cCntbit[0];
		iSrt	= 0;
		iEnd	= 1;
		iLen = DATA_NUMBER_LENGTH;
	}
	else {
		if (bFirst) {
			cEX  = tOut->cS;
			iSrt = tOut->iRan[idx];
			iEnd = tOut->iRan[idx+1];
			iLen = DATA_SQUARE_LENGTH;
		}
		else {
			cEX  = tOut->cF;
			iSrt = tOut->iRan2[idx];
			iEnd = tOut->iRan2[idx+1];
			iLen = DATA_NUMBER_LENGTH;
		}
	}

	// -E(x_i) 결과 저장 (i=1~l) l=DATA_SQUARE_LENGTH or DATA_NUMBER_LENGTH
	paillier_ciphertext_t _cEXi[iLen];
	for (int j=0 ; j<iLen ; j++)
		#ifdef _DEBUG_INIT_1
		mpz_init(_cEXi[j].c);
		#else
		mpz_init2(_cEXi[j].c, 2*GMP_N_SIZE);
		#endif

	#ifdef _DEBUG_SecureBitDecomposition
	// l = 2^(-1) mod N
	DebugOut("l=2^(-1) mod N\t\t\t", tPre->pL.m, idx);
	// N-1
	DebugOut("N\t\t\t\t", mPubKey->n, idx);
	DebugOut("N-1\t\t\t\t", tPre->pN1.m, idx);
	#endif

	for (int i=iSrt ; i<iEnd ; i++) {

		#ifdef _DEBUG_SecureBitDecomposition
		// parameter : E(x)
		//DebugOut("input E(x)\t\t\t", cEX[i].c, idx);
		DebugDec("input x\t\t\t", &cEX[i], idx);
		#endif

		if (tSkle->bSkle)
			cEXi = tSkle->cCntbit[0];
		else
			cEXi = bFirst ? tOut->cSbit[i] : tOut->cFbit[i];

		while (1) {
			// T = E(x)
			mpz_set(cT.c, cEX[i].c);
			#ifdef _DEBUG_SecureBitDecomposition
			//DebugOut("T (copied)", cT.c, idx);
			DebugDec("Decrypted T = E(x) (copied)\t", &cT, idx); // problem
			#endif

			for (int j=0 ; j<iLen ; j++) {
				#ifdef _DEBUG_SecureBitDecomposition
				printf("\n<<<  %03d - th bit (Bit-Decomposition) >>>\n\n", j);
				#endif

				// compute the j-th LSB cf E(x)
				//EncryptedLSB(&cEXi[j], &cT, tPre, idx, tRecv, tSend, ucSendPtr);
				EncryptedLSB(&cEXi[j], &cT, tPre, idx, tRecv);
				#ifdef _DEBUG_SecureBitDecomposition
				//DebugOut(j+"-th bit\t\t\t", cEXi[j].c, idx);
				DebugDec("j-th bit\t\t\t", &cEXi[j], idx);
				#endif

				// E(x_i)^(N-1) = E(-x_i)
				paillier_exp(mPubKey, &_cEXi[j], &cEXi[j], &(tPre->pN1));
				#ifdef _DEBUG_SecureBitDecomposition
				//DebugOut("E(-x_i)\t\t\t", _cEXi[j].c, idx);
				DebugDec("E(-x_i)\t\t\t", &_cEXi[j], idx);
				#endif

				// Z = T * E(-x_i)
				paillier_mul(mPubKey, &cZ, &cT, &_cEXi[j]);
				#ifdef _DEBUG_SecureBitDecomposition
				//DebugOut("Z = T * E(-x_i)\t\t", cZ.c, idx);
				DebugDec("Z = T * E(-x_i)\t\t", &cZ, idx);
				#endif

				// T = Z^l (right shift)
				paillier_exp(mPubKey, &cT, &cZ, &(tPre->pL));
				#ifdef _DEBUG_SecureBitDecomposition
				//DebugOut("T = Z^l (right shift)\t\t", T.c, idx);
				DebugDec("T = Z^l (right shift)\t\t", &cT, idx);
				#endif
			}

			#ifdef _DEBUG_SecureBitDecomposition
			// parameter : E(x), < E(x_0), E(x_1), ... , E(x_m-1) >
			DebugDec("x (Decrypted Parameter E(x))\t", &cEX[i], idx);
			printf("[DH] < x_m-1, ..., x_0, x_1 > (Decrypted Parameters)\n");
			for (int j=iLen-1 ; j>=0 ; j--)
				DebugDecBit(&cEXi[j]);
			printf("\n");
			#endif

			//iGamma = SVR(&cEX[i], cEXi, tPre, idx, tRecv, tSend, ucSendPtr);
			iGamma = SVR(&cEX[i], cEXi, iLen, tPre, idx, tRecv);
			if (iGamma==1) {
				if (tSkle->bSkle) {
					for (int j=0 ; j<iLen ; j++)
						// E(x_j) = 1 - E(x_j)
						paillier_mul(mPubKey, &(tSkle->cCntbit[1][j]), &(tPre->c1), &_cEXi[j]);

					#ifdef _DEBUG_SecureBitDecomposition
					// parameter : E(x), < E(x_0), E(x_1), ... , E(x_m-1) >
					DebugDec("x (Decrypted Parameter E(x))\t", &cEX[i], idx);
					printf("[DH] < x_m-1, ..., x_0, x_1 > (Decrypted Parameters)\n");
					for (int j=iLen-1 ; j>=0 ; j--)
						DebugDecBit(&(tSkle->cCntbit[1][j]));
					printf("\n");
					#endif
				}
				else {
					if (bFirst) {
						for (int j=0 ; j<iLen ; j++)
							// E(x_j) = 1 - E(x_j)
							paillier_mul(mPubKey, &cEXi[j], &(tPre->c1), &_cEXi[j]);

						#ifdef _DEBUG_SecureBitDecomposition
						// parameter : E(x), < E(x_0), E(x_1), ... , E(x_m-1) >
						DebugDec("x (Decrypted Parameter E(x))\t", &cEX[i], idx);
						printf("[DH] < x_m-1, ..., x_0, x_1 > (Decrypted Parameters)\n");
						for (int j=iLen-1 ; j>=0 ; j--)
							DebugDecBit(&cEXi[j]);
						printf("\n");
						#endif
					}
				}

				// checking the allocated size  of GMP
				#ifdef _DEBUG_Assert
				assert(cT.c->_mp_alloc == 2*GMP_N_SIZE);
				assert(cZ.c->_mp_alloc == 2*GMP_N_SIZE*2);
				assert(_cEXi[0].c->_mp_alloc == 2*GMP_N_SIZE);
				#endif

				break;
			}

			printf("\n<<<  ***  FAIL to SVR !!!  *** >>>\n\n");
		}
	}
}

//void PaillierCrypto::EncryptedLSB(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cT, pre_t* tPre, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::EncryptedLSB(paillier_ciphertext_t* cRes, paillier_ciphertext_t* cT, pre_t* tPre, unsigned short idx, thp_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bY[ENC_SIZE]={0,};
	short sLen;

	paillier_plaintext_t  pT;
	paillier_ciphertext_t cRY, c1;
	int iOdd;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pT.m, cRY.c, NULL);
	#else
	mpz_init2(pT.m, 2*GMP_N_SIZE);
	mpz_init2(cRY.c, 2*GMP_N_SIZE*2);
	#endif

	// parameter : Enrypted T
	#ifdef _DEBUG_EncryptedLSB
	//DebugOut("input T\t\t\t", cT->c, idx);
	DebugDec("input T\t\t\t", cT, idx);
	#endif

	// 0 <= r < N  (Pick random number)
	mpz_urandomm(pT.m, state, mPubKey->n); 				// pT = pR
	#ifdef _DEBUG_EncryptedLSB
	DebugOut("r (random number)\t\t", pT.m, idx);
	#endif

	// E(r)
	paillier_enc(&cRY, mPubKey, &pT, paillier_get_rand_devurandom); 	// 여기서부터는 cRY = cR
	#ifdef _DEBUG_EncryptedLSB
	//DebugOut("E(r)\t\t\t\t", cRY.c, idx);
	DebugDec("r\t\t\t\t", &cRY, idx);
	#endif

	// Y = T * E(r)
	paillier_mul(mPubKey, &cRY, cT, &cRY);				// 여기서부터는 cRY = cY
	#ifdef _DEBUG_EncryptedLSB
	//DebugOut("Y = T * E(r)\t\t\t", cRY.c, idx);
	DebugDec("y = x + r\t\t\t", &cRY, idx);
	#endif

	// send Y
	paillier_ciphertext_to_bytes(bY, ENC_SIZE, &cRY);
	SetSendMsg(ucSendPtr, bY, idx, COM_LSB, ENC_SIZE);
	#ifdef _DEBUG_EncryptedLSB
    DebugCom("Y = T * E(r) (Hex):\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

//    mSendMtx->lock();
//    mSendQueue->push(ucSendPtr);
//    mSendMtx->unlock();
//    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?

	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	// determine whether random r is odd
	iOdd = mpz_odd_p(pT.m);
	#ifdef _DEBUG_EncryptedLSB
	DebugOut("r is odd?\t\t\t", iOdd, idx);
	#endif

	// E(1) (if r is odd)
	if (iOdd) {
		#ifdef _DEBUG_INIT_1
		mpz_init(c1.c);
		#else
		mpz_init2(c1.c, 3*GMP_N_SIZE);
		#endif

		mpz_set_ui(pT.m, 1);						// 임시변수로서 pT를 사용함.
		paillier_enc(&c1, mPubKey, &pT, paillier_get_rand_devurandom);	// problem
		#ifdef _DEBUG_EncryptedLSB
		//DebugOut("E(1)", c1->c, idx);
		DebugDec("r\t\t\t\t", &c1, idx);
		#endif

		// N-1
		#ifdef _DEBUG_EncryptedLSB
		DebugOut("N\t\t\t\t", mPubKey->n, idx);	DebugOut("N-1\t\t\t\t", tPre->pN1.m, idx);
		#endif
	}

	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive alpha
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_EncryptedLSB
	DebugCom("alpha' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	#ifdef _DEBUG_Assert
    assert((*(ucRecvPtr+2)==COM_LSB)&&(sLen==ENC_SIZE));
	#endif

	paillier_ciphertext_from_bytes(cRes, ucRecvPtr+HED_SIZE, ENC_SIZE);
	#ifdef _DEBUG_EncryptedLSB
	//DebugOut("alpha' (Copied)", cAlpha.c, idx);
	DebugDec("alpha' (Copied)\t\t", cRes, idx);
	#endif

	// (r:even) E(x_i) = alpha
	// (r:odd)  E(x_i) = E(1) * alpha^(N-1)
	if (iOdd) {
		// alpha^(N-1)
		paillier_exp(mPubKey, &cRY, cRes, &(tPre->pN1));			// 임시 변수로서 cRY를 사용함.
		#ifdef _DEBUG_EncryptedLSB
		//DebugOut("alpha^(N-1)", cRY.c, idx);
		DebugDec("alpha^(N-1)\t\t\t", &cRY, idx);
		#endif

		// E(x_i) = E(1) * alpha^(N-1)
		paillier_mul(mPubKey, cRes, &c1, &cRY);
		#ifdef _DEBUG_EncryptedLSB
		//DebugOut("E(x_i) = E(1) * alpha^(N-1)", cRes->c, idx);
		DebugDec("E(x_i) = E(1) * alpha^(N-1)\t", cRes, idx);
		#endif
	}

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert( pT.m->_mp_alloc <= 2*GMP_N_SIZE);
	assert(cRY.c->_mp_alloc == 2*GMP_N_SIZE*2);
	if (iOdd)
		assert(c1.c->_mp_alloc <= 3*GMP_N_SIZE);
	#endif

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return;						// 2N*2
}

//int PaillierCrypto::SVR(paillier_ciphertext_t* cEX, paillier_ciphertext_t* cEXi, int len, pre_t* tPre, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
int PaillierCrypto::SVR(paillier_ciphertext_t* cEX, paillier_ciphertext_t* cEXi, int len, pre_t* tPre, unsigned short idx, thp_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bW[ENC_SIZE]={0,};
	short sLen;

	paillier_plaintext_t pT;
	paillier_ciphertext_t cUVW, cT;
	int iGamma;

	#ifdef _DEBUG_INIT_1
	mpz_inits(pT.m, cUVW.c, cT.c, NULL);
	#else
	mpz_init2(pT.m, GMP_N_SIZE+1);
	mpz_init2(cUVW.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE);
	#endif

	// parameter : E(x)
	#ifdef _DEBUG_SVR
	DebugOut("input E(x)\t\t\t", cEX->c, idx);
	DebugDec("input x\t\t\t", cEX, idx);
	#endif

	#ifdef _DEBUG_SVR
	// parameter : < E(x_0), E(x_1), ... , E(x_m-1) >
    printf("input < x_0, x_1, ..., x_m-1 >\n");
    for (int j=len-1 ; j>=0 ; j--) {
		DebugDecBit(&cEXi[j]);
		if (j%8==0)		printf("\t");
		if (j%64==0)	printf("\n");
	}
	#endif

	// U = SUM_i=0~m-1 (E(x_i))^(2^i)
    mpz_set_ui(pT.m, 1);
    mpz_set_ui(cUVW.c, 1);
	for (int i=0 ; i<len ; i++) {
		paillier_exp(mPubKey, &cT, &cEXi[i], &pT);				// 여기부터 pT = 2^i
		paillier_mul(mPubKey, &cUVW, &cUVW, &cT);				// 여기부터 cUVW = cU
		#ifdef _DEBUG_SVR
		//DebugOut("U = SUM_i=0~m-1 (E(x_i))^(2^i)", cUVW.c, idx);
		DebugDec("U = SUM_i=0~m-1 (E(x_i))^(2^i)", &cUVW, idx);
		#endif

		mpz_mul_2exp(pT.m, pT.m, 1);
	}

	// N-1
	#ifdef _DEBUG_SVR
	DebugOut("N\t\t\t\t", mPubKey->n, idx);	DebugOut("N-1\t\t\t\t", tPre->pN1.m, idx);
	#endif

	// V = U * E(x)^(N-1)
	paillier_exp(mPubKey, &cT, cEX, &(tPre->pN1));
	paillier_mul(mPubKey, &cUVW, &cUVW, &cT);				// 여기부터 cUVW = cV
	#ifdef _DEBUG_SVR
	//DebugOut("V = U * E(-x)", cUVW.c, idx);
	DebugDec("V = U * E(-x)\t\t\t", &cUVW, idx);
	#endif

	// 0 <= r < N  (Pick random number)
	mpz_urandomm(pT.m, state, mPubKey->n);					// 여기부터 pT = r
	#ifdef _DEBUG_SVR
	DebugOut("r (random number)\t\t", pT.m, idx);
	#endif

	// W = V^r
	paillier_exp(mPubKey, &cUVW, &cUVW, &pT);				// 여기부터 cUVW = cW
	#ifdef _DEBUG_SVR
	//DebugOut("W = V^r", cUVW.c, idx);
	DebugDec("W = V^r\t\t\t", &cUVW, idx);
	#endif

    // send W
	paillier_ciphertext_to_bytes(bW, ENC_SIZE, &cUVW);
	SetSendMsg(ucSendPtr, bW, idx, COM_SVR, ENC_SIZE);
	#ifdef _DEBUG_SVR
    DebugCom("W (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
	#endif

//    mSendMtx->lock();
//    mSendQueue->push(ucSendPtr);
//    mSendMtx->unlock();
//    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?

	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// receive ooo
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_SVR
	DebugCom("Gamma' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	// a' = E(a+r_a)
	#ifdef _DEBUG_Assert
    assert((*(ucRecvPtr+2)==COM_SVR)&&(sLen==1));
	#endif

    iGamma = *(ucRecvPtr+HED_SIZE);

	// E( (a+r_a)^2 ) = E( a^2 + 2*a*r_a + r_a^2 )
	#ifdef _DEBUG_SVR
	printf("[DH-%03d] Gamma \t\t\t\t : %d\n", idx, iGamma);
	#endif

	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(  pT.m->_mp_alloc ==   GMP_N_SIZE);
	assert(cUVW.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(  cT.c->_mp_alloc == 2*GMP_N_SIZE);
	#endif

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;

	return iGamma;
}

void PaillierCrypto::SkLE_s_Init(out_t* tOut, skle_t* tSkle, pre_t* tPre, bool bFirst, unsigned short idx)
{
	tSkle->bSkle = true;

	// Res'_i = E(0), Can'_i = E(1), tRes'_i for i=1~n
	int iNum = bFirst ? tOut->iRan[idx+1]-tOut->iRan[idx] : tOut->iRan2[idx+1]-tOut->iRan2[idx];
	for (int i=0 ; i<iNum ; i++) {
//		mpz_set_ui(tSkle->cIRes[i].c, 1);
//		mpz_set(tSkle->cCan[i].c, tPre->c1.c);
		mpz_set(tSkle->cC[i].c, tPre->c1.c);
		mpz_set_ui(tSkle->cK[i].c, 1);
	}

	#ifdef _DEBUG_SkLE_s
	printf("[DH-%03d] Can \t\t\t : ", idx);
	for (int i=0 ; i<iNum ; i++)
		DebugDecBit(&(tSkle->cC[i]));
	printf("\n");

	printf("[DH-%03d] Res \t\t\t : ", idx);
	for (int i=0 ; i<iNum ; i++)
		DebugDecBit(&(tSkle->cK[i]));
	printf("\n");
	#endif

	return;
}

//void PaillierCrypto::SkLE_s_Main1(out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::SkLE_s_Main(out_t* tOut, skle_t* tSkle, pre_t* tPre, bool bFirst, unsigned short idx, thp_t* tRecv)
{
	int iNum = bFirst ? tOut->iRan[idx+1]-tOut->iRan[idx] : tOut->iRan2[idx+1]-tOut->iRan2[idx];
	int iLen = bFirst ? DATA_SQUARE_LENGTH : DATA_NUMBER_LENGTH;

	paillier_ciphertext_t  *cE, *cC, *cK;
	paillier_ciphertext_t  cS, cP;
	paillier_ciphertext_t  cU[iNum], cT;
	paillier_ciphertext_t  cM, cD, cA, cB, cG;

	#ifdef _DEBUG_INIT_1
	mpz_inits(cS.c, cP.c, cT.c, cM.c, cD.c, cA.c, cB.c, cG.c, NULL);
	for (int i=0 ; i<iNum ; i++)
		mpz_inits(cU[i].c, NULL);
	#else
	mpz_init2(cS.c, 2*GMP_N_SIZE*2);
	mpz_init2(cP.c, 2*GMP_N_SIZE*2);
	mpz_init2(cT.c, 2*GMP_N_SIZE*2);
	mpz_init2(cM.c, 2*GMP_N_SIZE*2);
	mpz_init2(cD.c, 2*GMP_N_SIZE*2);
	mpz_init2(cA.c, 2*GMP_N_SIZE*2);
	mpz_init2(cB.c, 2*GMP_N_SIZE*2);
	mpz_init2(cG.c, 2*GMP_N_SIZE*2);
	for (int i=0 ; i<iNum ; i++) {
		mpz_init2(cU[i].c, 2*GMP_N_SIZE*2);
	}
	#endif

	int SrtDat, EndDat;
	int iTidx ;

	SrtDat = bFirst ? tOut->iRan[idx]  : tOut->iRan2[idx];			// 각 thread의 data 시작점
	EndDat = SrtDat + iNum;												// 각 thread의 data 종료점
//	cC = tSkle->cC;
//	cK = tSkle->cK;




	for (int j=iLen-1 ; j>=0 ; j++) {			// j = l-1 ~ 0


		///////////////////////// Step 1 /////////////////////////

		// initialize S
		mpz_set_ui(cS.c, 1);

		for (int i=SrtDat ; i<EndDat ; i++) {			// i = 각 thread의 data 시작점 ~ 종료점
			// iTidx = bFirst ? tOut->iRan[idx]+i  : tOut->iRan2[idx]+i;
			cE = bFirst ? &(tOut->cSbit[i][j]) : &(tOut->cFbit[i][j]);
			cC = &(tSkle->cC[i]);
			cK = &(tSkle->cK[i]);

			#ifdef _DEBUG_SkLE_s
			printf("[DH-%03d] *** (step 1) [%d, %d] P ***\n", (int)idx, i, j);
			DebugDec("E[i,j]\t\t\t", cE, idx);
			DebugDec("C[i]\t\t\t\t", cC, idx);
			DebugDec("K[i]\t\t\t\t", cK, idx);
			#endif

			// U[i] = SM(C[i], e[i,j])
			//SecMul(&cTmp1i, &(tSkle->cCan[i]), &cE[j], idx, tRecv);
			SecMul(cU+i-SrtDat, cC, cE, idx, tRecv);
			#ifdef _DEBUG_SkLE_s
			DebugDec("U[i] = SM(C[i], e[i,j]) \t", cU+i-SrtDat, idx);
			#endif

			// P = K[i] * U[i]
			paillier_mul(mPubKey, &cP, cU+i-SrtDat, cK);
			#ifdef _DEBUG_SkLE_s
			DebugDec("P = K[i] * T[i] \t\t", &cP, idx);
			#endif

			// S = Sum_(i=1~n) P[i]
			//paillier_mul(mPubKey, &(tOut->cTCnt[idx]), &(tOut->cTCnt[idx]), &(tSkle->cTRes[i]));
			paillier_mul(mPubKey, &cS, &cS, &cP);
			#ifdef _DEBUG_SkLE_s
			DebugDec("S = Sum_(i=1~n) P[i] \t", &cS, idx);
			#endif
		}

		///////////////////////// Step 2 /////////////////////////

//		// thread 일시정지
//		ulSync.lock();
//		if (iNum%2 == 0) 		{	pSyncCnt2 = &(tSync->c2);	pSyncCnt1 = &(tSync->c1);	}
//		else 					{	pSyncCnt2 = &(tSync->c1);	pSyncCnt1 = &(tSync->c2);	}
//
//		(*pSyncCnt2)++;
//		if (*pSyncCnt2 < NUM_MAIN_THREAD) {
//			tSync->cv.wait(ulSync, [&] {return *pSyncCnt2>=NUM_MAIN_THREAD;});
//		}
//		else {
//			//oCrypto->SkLE_s_Cmp(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv, tSend, ucSendBuf);
//			oCrypto->SkLE_s_Cmp(tOut, &tSkle, &tPre, iBit-1, bFirst, idx, tRecv);
//			tSync->cv.notify_all();		*pSyncCnt1=0;
//		}
//		ulSync.unlock();


		// 각 thread의 S 합산


		///////////////////////// Step 3 /////////////////////////

		// SBD of S
		SecBitDecomp(tOut, tSkle, tPre, bFirst, idx, tRecv);

		// SCI

		// temp code: E(M)=E(1),E(D)=E(1)
		mpz_set(cM.c, tPre->c1.c);
		mpz_set(cD.c, tPre->c1.c);



		///////////////////////// Step 4 /////////////////////////

		// Alpha, Beta, Gamma 계산

		// Alpha = SM(D, M)
		//SecMul(&cTmp1i, &(tSkle->cCan[i]), &cE[j], idx, tRecv);
		SecMul(&cA, &cD, &cM, idx, tRecv);
		#ifdef _DEBUG_SkLE_s
		DebugDec("Alpha = SM(D, M) \t\t", &cA, idx);
		#endif

		// Beta = E(1) * E(D)^(N-1) * E(Alpha)
		paillier_exp(mPubKey, &cT, &cD, &(tPre->pN1));
		paillier_mul(mPubKey, &cT, &(tPre->c1), &cT);
		paillier_mul(mPubKey, &cB, &cT, &cA);
		#ifdef _DEBUG_SkLE_s
		DebugDec("Beta = E(1)*E(D)^(N-1)*E(A)", &cB, idx);
		#endif

		// Gamma = E(D) * E(Alpha)^(N-2)
		paillier_exp(mPubKey, &cT, &cA, &(tPre->pN2));
		paillier_mul(mPubKey, &cG, &cD, &cT);
		#ifdef _DEBUG_SkLE_s
		DebugDec("Gamma = E(D) * E(A)^(N-2) ", &cG, idx);
		#endif

		// thread 재시작



		for (int i=0 ; i<iNum ; i++) {
			//iTidx = bFirst ? tOut->iRan[idx]+i  : tOut->iRan2[idx]+i;
			cE = bFirst ? &(tOut->cSbit[i][j]) : &(tOut->cFbit[i][j]);
			cC = &(tSkle->cC[i]);
			cK = &(tSkle->cK[i]);

			#ifdef _DEBUG_SkLE_s
			printf("[DH-%03d] *** (step 4) [%d, %d] C, K ***\n", (int)idx, i, j);
			#endif

			// T = SM(U[i], Beta)
			//SecMul(&cTmp1i, &(tSkle->cCan[i]), &cE[j], idx, tRecv);
			SecMul(&cT, cU+i-SrtDat, &cB, idx, tRecv);
			#ifdef _DEBUG_SkLE_s
			DebugDec("T = SM(U[i], Beta) \t", &cT, idx);
			#endif

			// K[i] = K[i] * T
			paillier_mul(mPubKey, cK, cK, &cT);
			#ifdef _DEBUG_SkLE_s
			DebugDec("K[i] = K[i] * T \t\t", cK, idx);
			#endif

			// T = SM(e[i,j], Gamma)
			//SecMul(&cTmp1i, &(tSkle->cCan[i]), &cE[j], idx, tRecv);
			SecMul(&cT, cE, &cG, idx, tRecv);
			#ifdef _DEBUG_SkLE_s
			DebugDec("T = SM(e[i,j], Gamma) \t", &cT, idx);
			#endif

			// T = M * T
			paillier_mul(mPubKey, &cT, &cM, &cT);
			#ifdef _DEBUG_SkLE_s
			DebugDec("T = M * T \t\t\t", &cT, idx);
			#endif

			// C[i] = SM(C[i], T)
			//SecMul(&cTmp1i, &(tSkle->cCan[i]), &cE[j], idx, tRecv);
			SecMul(cC+i, cC, &cT, idx, tRecv);
			#ifdef _DEBUG_SkLE_s
			DebugDec("C[i] = SM(C[i], T) \t", cC, idx);
			#endif

			//iTidx++;
		}

		// print Res, Can for debugging
		//printf("[DH-%03d] %d-th bit Can : ", idx, j);	for (int i=0 ; i<iNum ; i++)	DebugDecBit(&(tSkle->cCan[i]));		printf("\n");
		//printf("[DH-%03d] %d-th bit Res : ", idx, j);	for (int i=0 ; i<iNum ; i++)	DebugDecBit(&(tSkle->cIRes[i]));		printf("\n");

	}


	for (int i=0 ; i<iNum ; i++) {			// i = 1 ~ n
		// K[i] = K[i] * C[i]
		paillier_mul(mPubKey, cK+i, cK+i, cC+i);
		#ifdef _DEBUG_SkLE_s
		DebugDec("K[i] = K[i] * C[i] \t", cK+i, idx);
		#endif
	}


	#ifdef _DEBUG_Assert
	// checking the allocated size  of GMP
	assert(cTmp1i.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(cTmp2i.c->_mp_alloc <= 2*GMP_N_SIZE*2);
	assert(cTmp3i.c->_mp_alloc <= 2*GMP_N_SIZE*2);
	#endif

	return;
}

//void PaillierCrypto::SkLE_s_Cmp(out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::SkLE_s_Cmp(out_t* tOut, skle_t* tSkle, pre_t* tPre, int bit, bool bFirst, unsigned short idx, thp_t* tRecv)
{
	int iNumThread, iCmpRes, iTmp ;

	iNumThread = bFirst ? NUM_MAIN_THREAD : NUM_MAIN_THREAD2;

	// initialize cCnt
	mpz_set_ui(tOut->cCnt.c, 1);
	// Cnt' = Sum_(i=1~n) tRes'_i
	for (int i=0 ; i<iNumThread ; i++)
		paillier_mul(mPubKey, &(tOut->cCnt), &(tOut->cCnt), &(tOut->cTCnt[i]));
	#ifdef _DEBUG_SkLE_s
	//DebugOut("Cnt' = Sum_(i=1~n) tRes'_i", tOut->cCnt.c, idx);
	DebugDec("Cnt' = Sum_(i=1~n) tRes'_i\t", &(tOut->cCnt), idx);
	#else
	//printf("%d \n", bit);
	printf("%d bit - ", bit);
	DebugDec("Cnt'\t\t", &(tOut->cCnt), idx);
	#endif

	// compute Secure Bit-Decomposition of Cnt
	//SecBitDecomp(tOut, &tSkle, &tPre, bFirst, idx, tRecv, tSend, ucSendBuf);
	SecBitDecomp(tOut, tSkle, tPre, bFirst, idx, tRecv);

	if (bFirst) {		// param_k=PARAM_Kbit
		iTmp = PARAM_K;
		for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
			tSkle->iKbit[i] = iTmp%2;
			iTmp >>= 1;
		}
		iTmp = PARAM_K;
	}
	else {				// param_k=1
		std::fill_n(tSkle->iKbit, DATA_NUMBER_LENGTH, 0);
		tSkle->iKbit[0] = 1;
		iTmp = 1;
	}

//	iCmpRes = Comparison(&(tOut->cCmp), &(tOut->cCnt), tSkle->cCntbit[0], tSkle->cCntbit[1], iTmp, tSkle->iKbit, tPre, idx, tRecv);
//	tOut->iCmp = Comparison(&(tOut->cCmp), &(tOut->cCnt), tSkle->cCntbit[0], tSkle->cCntbit[1], iTmp, tSkle->iKbit, tPre, idx, tRecv);
	#ifdef _DEBUG_SkLE_s
//	DebugOut("Comparison result\t\t", iCmpRes, idx);
	DebugOut("Comparison result\t\t", tOut->iCmp, idx);
	//DebugOut("Comparison data\t\t", cCmp->c, idx);
	DebugDec("Comparison data\t\t", &(tOut->cCmp), idx);
	#endif

	if ((iCmpRes==0) || (bit<0)) {
		if (iCmpRes == 0)						tOut->iCmp = 1;		// cTRes;
		if (bit < 0) 							tOut->iCmp = 2;		// cIRes;
	}
	else										tOut->iCmp = 0;

	return;
}

//int PaillierCrypto::Comparison(paillier_ciphertext_t* cGamma, paillier_ciphertext_t* cCnt, paillier_ciphertext_t* cCntbit,
//		paillier_ciphertext_t* cCntCom, int iK, int* iKbit, pre_t *tPre, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)

//int PaillierCrypto::Comparison(paillier_ciphertext_t* cM, paillier_ciphertext_t* cD, paillier_ciphertext_t* cCnt, paillier_ciphertext_t* cCntbit,
//		paillier_ciphertext_t* cCntCom, int iK, int* iKbit, pre_t *tPre, unsigned short idx, thp_t* tRecv)
//{
//	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
//	unsigned char* ucRecvPtr=NULL;
//	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
//	unsigned char bX[ENC_SIZE]={0,};
//	short sLen;
//
//	paillier_plaintext_t p0, pR, pT ;
//	paillier_ciphertext_t cW[DATA_NUMBER_LENGTH], cX[DATA_NUMBER_LENGTH], cY[DATA_NUMBER_LENGTH] ;
//	paillier_ciphertext_t cGamma1, cBeta, cP, cT, cM, cD, c1 ;
//	int iAlpha, iRes;
//	unsigned int uiPhi[DATA_SIZE], uiTmp, k;
//
//	#ifdef _DEBUG_INIT_1
//	mpz_inits(p0.m, pR.m, pT.m, cGamma1.c, cBeta.c, cP.c, cT.c, cM->c, cD->c, c1.c, NULL);
//	for (int j=0 ; j<DATA_NUMBER_LENGTH ; j++)
//		mpz_inits(cW[j].c, cX[j].c, cY[j].c, NULL);
////		mpz_init(cU[j].c);
//	#else
//	mpz_init2(p0.m, GMP_N_SIZE+1);
//	mpz_init2(pR.m, GMP_N_SIZE+1);
//	mpz_init2(pT.m, GMP_N_SIZE+1);
//	mpz_init2(cGamma.c, 3*GMP_N_SIZE);
//	mpz_init2(cP.c, 3*GMP_N_SIZE);
//	mpz_init2(cT.c, 3*GMP_N_SIZE);
//	mpz_init2(cM.c, 3*GMP_N_SIZE);
//	mpz_init2(cD.c, 3*GMP_N_SIZE);
//	mpz_init2(cV.c, 3*GMP_N_SIZE);
//	mpz_init2(cB.c, 3*GMP_N_SIZE);
//	mpz_init2(c1.c, 3*GMP_N_SIZE);
//	for (int j=0 ; j<DATA_NUMBER_LENGTH ; j++)
//		mpz_init2(cU[j].c, 3*GMP_N_SIZE);
//	#endif
//
//	srand((unsigned int)time(0));
//
//	#ifdef _DEBUG_Comparison
//	DebugDec("input Cnt'\t\t\t", cCnt, idx);
//	DebugOut("input Parameter k\t\t", iK, idx);
//	printf("[DH] < Cnt_m-1, ..., Cnt_1, Cnt_0 > (Decrypted Count)\n");
//	for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {
//		DebugDecBit(&cCntbit[j]);
//		if (j%8==0)		printf("\t");
//		if (j%64==0)	printf("\n");
//	}
//	printf("\n");
//	printf("[DH] < ~Cnt_m-1, ..., ~Cnt_1, ~Cnt_0 > (Decrypted ~Count)\n");
//	for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {
//		DebugDecBit(&cCntCom[j]);
//		if (j%8==0)		printf("\t");
//		if (j%64==0)	printf("\n");
//	}
//	printf("\n");
//	printf("[DH] < k_m-1, ..., k_1, k_0 > (Decrypted Parameter k)\n");
//	for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {
//		printf("%d", iKbit[j]);
//		if (j%8==0)		printf("\t");
//		if (j%64==0)	printf("\n");
//	}
//	printf("\n");
//	#endif
//
//	mpz_set_ui(p0.m, 0);
//	mpz_set_ui(cP.c, 1);
//	mpz_set_ui(cY[DATA_NUMBER_LENGTH].c, 1);			//?? y_l'=E(0) 맞나?
//
//	mpz_set_ui(pT.m, 1);						// 임시변수로서 pT를 사용함.
//	paillier_enc(&c1, mPubKey, &pT, paillier_get_rand_devurandom);	// problem
//	#ifdef _DEBUG_EncryptedLSB
//	//DebugOut("E(1)", c1->c, idx);
//	DebugDec("r\t\t\t\t", &c1, idx);
//	#endif
//
//	// alpha = {0, 1}
//	iAlpha = rand() % 2;
//	#ifdef _DEBUG_Comparison
//	DebugOut("alpha\t\t\t\t\t\t\t", iAlpha, idx);
//	#endif
//
//	for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {			//?? 첫번째와 두번째 길이가 달라질 것 같은데.
//		if (iAlpha==0) {
//			// E(T) = E(~Cnt_j * k_j) : E(0) or E(~Cnt_j)
//			if (iKbit[j]==0)		paillier_enc(&cT, mPubKey, &p0, paillier_get_rand_devurandom);
//			else					mpz_set(cT.c, cCntCom[j].c);
//			#ifdef _DEBUG_Comparison
//			DebugOut("Parameter k_j\t\t\t\t\t\t", iKbit[j], idx);
//			DebugDec("E(T) = E(~Cnt_j * k_j) : E(0) or E(~Cnt_j) (Alpha=0)\t", &cT, idx);
//			#endif
//		}
//		else {
//			// E(T) = E(Cnt_j * ~k_j) : E(Cnt_j) or E(0)
//			if (iKbit[j]==0)		mpz_set(cT.c, cCntbit[j].c);
//			else					paillier_enc(&cT, mPubKey, &p0, paillier_get_rand_devurandom);
//			#ifdef _DEBUG_Comparison
//			DebugOut("Parameter k_j\t\t\t\t\t\t", iKbit[j], idx);
//			DebugDec("E(T) = E(Cnt_j * ~k_j) : E(Cnt_j) or E(0) (Alpha=1)\t", &cT, idx);
//			#endif
//		}
//
//		// E(W) = E(W)^r
//		mpz_urandomm(pR.m, state, mPubKey->n);
//		paillier_exp(mPubKey, &cW[j], &cT, &pR);
//		#ifdef _DEBUG_Comparison
//		DebugDec("E(W_j) = E(T)^r \t\t\t\t", &cW[j], idx);
//		#endif
//
//		// E(X_j) = E(S_j xor k_j) : E(Cnt_j) or E(~Cnt_j)
//		if (iKbit[j]==0)		mpz_set(cX[j].c, cCntbit[j].c);
//		else					mpz_set(cX[j].c, cCntCom[j].c);
//		#ifdef _DEBUG_Comparison
//		DebugOut("Parameter k_j\t\t\t\t\t\t", iKbit[j], idx);
//		DebugDec("E(X_j) = E(S_j xor k_j) : E(Cnt_j) or E(~Cnt_j)\t", &cX[j], idx);
//		#endif
//
//		// E(Y_j) = E(Y_j+1)^r * E(X_j)
//		mpz_urandomm(pR.m, state, mPubKey->n);
//		paillier_exp(mPubKey, &cY[j], &cY[j+1], &pR);
//		paillier_mul(mPubKey, &cY[j], &cY[j], &cX[j]);
//		#ifdef _DEBUG_Comparison
//		DebugDec("E(U_j) = E(U_j+1)^r * E(T)\t\t\t\t", &cY[j], idx);
//		#endif
//	}
//
//	// E(Gamma) = SZP( E(x_all) )
//	mpz_set_ui(cGamma1.c, 1);			//?? E(Gamma1)=E(0) 맞나?
//	#ifdef _DEBUG_Comparison
//	DebugDec("E(Gamma) = SZP( E(x_all) \t\t\t\t\t", &cGamma1, idx);
//	#endif
//
//	for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {
//		// E(T) = E(y_j) * E(N-1)
//		paillier_mul(mPubKey, &cT, &cY[j], &(tPre->cN1));
//		#ifdef _DEBUG_Comparison
//		DebugDec("E(T) = E(Y_j) * E(N-1)\t\t\t\t", &cT, idx);
//		#endif
//
//		// E(X_j) = E(T)^r_j * E(W_j) : 예전 메모리 활용
//		mpz_urandomm(pR.m, state, mPubKey->n);
//		paillier_exp(mPubKey, &cT, &cT, &pR);
//		paillier_mul(mPubKey, &cX[j], &cT, &cW[j]);
//		#ifdef _DEBUG_Comparison
//		DebugDec("E(X_j) = E(T)^r_j * E(W_j)\t\t\t\t", &cX[j], idx);
//		#endif
//	}
//
//	//	// E(V) = E(r(Cnt-k))
//	//	if (iK!=1)	paillier_mul(mPubKey, &cV, cCnt, &(tPre->cNk));		// data의 SkLE_s
//	//	else 		paillier_mul(mPubKey, &cV, cCnt, &(tPre->cN1));		// frequency의 SkLE_s (k=1)
//	//	#ifdef _DEBUG_Comparison
//	//	DebugDec("E(Cnt-k)\t\t\t\t\t\t", &cV, idx);
//	//	#endif
//	//	mpz_urandomm(pR.m, state, mPubKey->n);
//	//	paillier_exp(mPubKey, &cV, &cV, &pR);
//	//	#ifdef _DEBUG_Comparison
//	//	DebugDec("E(V) = E(r(Cnt-k))\t\t\t\t\t", &cV, idx);
//	//	#endif
//
//	// ** generate random permutation (Fisher?Yates shuffle Algorithm) ** //
//
//	// Phi permutation
//	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++)
//		uiPhi[i] = i;
//	for (int i=DATA_NUMBER_LENGTH-1 ; i>0 ; i--) {
//		k = rand() % (i+1);
//		uiTmp = uiPhi[i];
//		uiPhi[i] = uiPhi[k];
//		uiPhi[k] = uiTmp;
//	}
//
//	#ifdef _DEBUG_Comparison
//	printf("[DH-%03d] Phi permutation : ", idx);
//	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++)
//		printf(" %d ", uiPhi[i]);
//	std::cout << std::endl;
//	#endif
//
//	for (int j=DATA_NUMBER_LENGTH-1 ; j>=0 ; j--) {
//	    // send Phi( X_all )
//		paillier_ciphertext_to_bytes(bX, ENC_SIZE, &cX[uiPhi[j]]);
//		SetSendMsg(ucSendPtr, bX, idx, COM_SVR, ENC_SIZE);
//		#ifdef _DEBUG_SVR
//	    DebugCom("X (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
//		#endif
//
////    mSendMtx->lock();
////    mSendQueue->push(ucSendPtr);
////    mSendMtx->unlock();
////    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?
//
//		try {
//			mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
//		}
//		catch ( SocketException& e ) {
//			std::cout << "Exception: " << e.description() << std::endl;
//		}
//	}
//
//	// E(D) = E(1) * gamma^(N-1)
//    // 데이터 전송 직후 여유시간에 계산
//
//	// gamma^(N-1)
//	paillier_exp(mPubKey, &cD, &cG, &(tPre->pN1));			// 임시 변수로서 cRY를 사용함.
//	#ifdef _DEBUG_EncryptedLSB
//	//DebugOut("Gamma^(N-1)", cD.c, idx);
//	DebugDec("Gamma^(N-1)\t\t\t", &cD, idx);
//	#endif
//
//	// E(D) = E(1) * gamma^(N-1)
//	paillier_mul(mPubKey, &cD, &c1, &cD);
//	#ifdef _DEBUG_EncryptedLSB
//	//DebugOut("E(D) = E(1) * gamma^(N-1)", cD->c, idx);
//	DebugDec("E(D) = E(1) * gamma^(N-1)\t", cD, idx);
//	#endif
//
//
//	ulRecv.lock();
//	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
//	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
//	ulRecv.unlock();
//
//	// receive E(beta)
//	sLen = Byte2Short(ucRecvPtr+HED_LEN);
//	#ifdef _DEBUG_Comparison
//    DebugCom("Beta' (Hex)\t\t\t\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
//	#endif
//
//	#ifdef _DEBUG_Assert
//    assert((*(ucRecvPtr+2)==COM_CMP1)&&(sLen==1));
//	#endif
//
//	paillier_ciphertext_from_bytes(&cBeta, ucRecvPtr+HED_SIZE, ENC_SIZE);
//	#ifdef _DEBUG_EncryptedLSB
//	//DebugOut("beta' (Copied)", cBeta.c, idx);
//	DebugDec("Beta' (Copied)\t\t", cBeta, idx);
//	#endif
//
//	if (iAlpha==0) {
//		mpz_set(cM->c, cBeta.c);
//		#ifdef _DEBUG_Comparison
//		DebugDec("E(M) = E(Beta) \t\t\t\t", cM, idx);
//		#endif
//	}
//	else {
//		// beta^(N-1)
//		paillier_exp(mPubKey, &cBeta, &cBeta, &(tPre->pN1));			// 임시 변수로서 cRY를 사용함.
//		#ifdef _DEBUG_EncryptedLSB
//		//DebugOut("Beta^(N-1)", cBeta.c, idx);
//		DebugDec("Beta^(N-1)\t\t\t", &cBeta, idx);
//		#endif
//
//		// E(M) = E(1) * beta^(N-1)
//		paillier_mul(mPubKey, &cBeta, &c1, &cBeta);
//		#ifdef _DEBUG_EncryptedLSB
//		//DebugOut("E(M) = E(1) * Beta^(N-1)", cBeta->c, idx);
//		DebugDec("E(M) = E(1) * beta^(N-1)\t", cBeta, idx);
//		#endif
//	}
//
//	// E(D) = E(1) * gamma^(N-1)
//    // (위에서 계산) 데이터 전송 직후 계산
//
//	// checking the allocated size  of GMP
//	#ifdef _DEBUG_Assert
//	assert(   p0.m->_mp_alloc ==   1);
//	assert(   pR.m->_mp_alloc ==   GMP_N_SIZE);
//	assert(   cP.c->_mp_alloc == 2*GMP_N_SIZE);
//	assert(   cT.c->_mp_alloc == 2*GMP_N_SIZE+1);
////	assert(   cV.c->_mp_alloc == 2*GMP_N_SIZE*2);
////	assert(cU[0].c->_mp_alloc == 2*GMP_N_SIZE*2);
//	#endif
//
//	// COM_CMP1 명령에서 (Cnt==k)일 경우에 대한 Recv Buff 초기화
//	// for loop의 COM_CMP2 명령 완료 후 break 되었을 경우에 대한 Recv Buff 초기화
//    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
//	ucRecvPtr[1] = 0xff;
//
//	return iRes;
//
//
//
////	if (iAlpha == 0) {
////		mpz_set(cGamma->c, cBeta.c);
////		#ifdef _DEBUG_Comparison
////		//DebugOut("(alpha=0) Gamma = E(Beta)", cGamma->c, idx);
////		DebugDec("(alpha=0) Gamma = E(Beta)\t\t\t\t", cGamma, idx);
////		#endif
////	}
////	else {		// if (iAlpha == 1)
////		paillier_exp(mPubKey, &cT, &cBeta, &(tPre->pN1));		// 임시변수로서 cTmp를 사용함.
////		paillier_mul(mPubKey, cGamma, &(tPre->c1), &cT);
////		#ifdef _DEBUG_Comparison
////		//DebugOut("(alpha=1) Gamma = E(1-Beta)", cGamma->c, idx);
////		DebugDec("(alpha=1) Gamma = E(1-Beta)\t\t\t\t", cGamma, idx);
////		#endif
////	}
////
////	#ifdef _DEBUG_Assert
////	assert(cBeta.c->_mp_alloc == 2*GMP_N_SIZE);
////	#endif
////
////
////
////
////	for (int i=0 ; i<DATA_NUMBER_LENGTH ; i++) {
////		// COM_CMP1 명령 완료에 대한 Recv Buff 초기화
////		// for loop의 COM_CMP2 명령 완료에 대한 Recv Buff 초기화 (break 되었을 경우는 제외)
////		memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
////		ucRecvPtr[1] = 0xff;
////
////		k = uiPhi[i];
////		#ifdef _DEBUG_Comparison
////		DebugOut("x-th bit\t\t\t\t\t\t", k, idx);
////		#endif
////
////		// send U'
////		paillier_ciphertext_to_bytes(bV, ENC_SIZE, &cU[k]);
////		SetSendMsg(ucSendPtr, bV, idx, COM_CMP2, ENC_SIZE);
////		#ifdef _DEBUG_Comparison
////		DebugCom("U_k' (Hex)\t\t\t\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
////		#endif
////
////		//    mSendMtx->lock();
////		//    mSendQueue->push(ucSendPtr);
////		//    mSendMtx->unlock();
////		//    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?
////
////		try {
////			mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
////		}
////		catch ( SocketException& e ) {
////			std::cout << "Exception: " << e.description() << std::endl;
////		}
////
////		ulRecv.lock();
////		tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
////		ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
////		ulRecv.unlock();
////
////		// beta'
////		sLen = Byte2Short(ucRecvPtr+HED_LEN);
////		#ifdef _DEBUG_Comparison
////		DebugCom("Beta' (Hex)\t\t\t\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
////		#endif
////
////		#ifdef _DEBUG_Assert
////		assert((*(ucRecvPtr+2)==COM_CMP2)&&((sLen==1)||(sLen==ENC_SIZE)));
////		#endif
////
////		if (sLen == ENC_SIZE) {			// Beta = E(0) or E(1)
////			paillier_ciphertext_t cBeta;
////			#ifdef _DEBUG_INIT_1
////			mpz_init(cBeta.c);
////			#else
////			mpz_init2(cBeta.c, 2*GMP_N_SIZE);
////			#endif
////			paillier_ciphertext_from_bytes(&cBeta, ucRecvPtr+HED_SIZE, ENC_SIZE);
////			#ifdef _DEBUG_Comparison
////			DebugOut("Beta' (Copied)\t\t\t\t\t\t", cBeta.c, idx);
////			DebugDec("Beta' (Copied)\t\t\t\t\t\t", &cBeta, idx);
////			#endif
////			iRes = 1;
////
////			if (iAlpha == 0) {
////				mpz_set(cGamma->c, cBeta.c);
////				#ifdef _DEBUG_Comparison
////				//DebugOut("(alpha=0) Gamma = E(Beta)", cGamma->c, idx);
////				DebugDec("(alpha=0) Gamma = E(Beta)\t\t\t\t", cGamma, idx);
////				#endif
////			}
////			else {		// if (iAlpha == 1)
////				paillier_exp(mPubKey, &cT, &cBeta, &(tPre->pN1));		// 임시변수로서 cTmp를 사용함.
////				paillier_mul(mPubKey, cGamma, &(tPre->c1), &cT);
////				#ifdef _DEBUG_Comparison
////				//DebugOut("(alpha=1) Gamma = E(1-Beta)", cGamma->c, idx);
////				DebugDec("(alpha=1) Gamma = E(1-Beta)\t\t\t\t", cGamma, idx);
////				#endif
////			}
////
////			#ifdef _DEBUG_Assert
////			assert(cBeta.c->_mp_alloc == 2*GMP_N_SIZE);
////			#endif
////
////			break;
////		}
////		else if ((sLen==1)&&(*(ucRecvPtr+HED_SIZE)==1)) {	// continue
////			iRes = 1;
////		}
////		else
////			assert(0);
////	}
//}

//void PaillierCrypto::CFTKD(out_t* tOut, in_t* tIn, pre_t* tPre, unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
void PaillierCrypto::CFTKD(out_t* tOut, in_t* tIn, pre_t* tPre, unsigned short idx, thp_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+ENC_SIZE]={0,};
	unsigned char bU[ENC_SIZE]={0,};
	short sLen;

	paillier_plaintext_t pR;
	paillier_ciphertext_t cTmp1, cW, cT[NUM_THREAD_DATA][CLASS_SIZE];
	unsigned int uiSigma[DATA_SIZE], uiPhi[DATA_SIZE][CLASS_SIZE];
	unsigned int uiTmp, k, l;
	int iTidx, iNum;

	#ifdef _DEBUG_INIT_1
//	mpz_inits(pTmp.m, pR.m, cTmp1.c, cW.c, NULL);
	mpz_inits(pR.m, cTmp1.c, cW.c, NULL);
	#else
//	mpz_init2(pTmp.m, 1);
	mpz_init2(pR.m, GMP_N_SIZE);
	mpz_init2(cTmp1.c, 2*GMP_N_SIZE*2);
	mpz_init2(cW.c, 2*GMP_N_SIZE);
	#endif
	for (int i=0 ; i<NUM_THREAD_DATA ; i++)
		for (int j=0 ; j<CLASS_SIZE ; j++)
			#ifdef _DEBUG_INIT_1
			mpz_init(cT[i][j].c);
			#else
			mpz_init2(cT[i][j].c, 2*GMP_N_SIZE*2);
			#endif

	srand((unsigned int)time(0));

	iNum = tOut->iRan[idx+1] - tOut->iRan[idx];

	// ** initialization ** //

	// compute (-c_j)' = E(N-c_j) for j=1~v
	#ifdef _DEBUG_CFTKD
	for (int j=0 ; j<CLASS_SIZE ; j++) {
		//DebugOut("(-c_j)' = E(N-c_j)", cCls[j].c, idx);
		DebugDec("(-c_j)' = E(N-c_j)", &(tPre->cCls[j]), idx);
	}
	#endif

	// ** generate random permutation (Fisher?Yates shuffle Algorithm) ** //

	// Sigma permutation (row)
	for (int i=0 ; i<iNum ; i++)
		uiSigma[i] = i;
	for (int i=iNum-1 ; i>0 ; i--) {
		k = rand() % (i+1);
		uiTmp = uiSigma[i];
		uiSigma[i] = uiSigma[k];
		uiSigma[k] = uiTmp;
	}

	#ifdef _DEBUG_CFTKD
	printf("[DH-%03d] Sigma permutation (row) : ", idx);
	for (int i=0 ; i<iNum ; i++)
		printf(" %d ", uiSigma[i]);
	std::cout << std::endl;
	#endif

	// Phi permutation (each column)
	for (int j=0 ; j<CLASS_SIZE ; j++)
		for (int i=0 ; i<iNum ; i++)
			uiPhi[i][j] = j;
	for (int i=0 ; i<iNum ; i++) {
		for (int j=CLASS_SIZE-1 ; j>0 ; j--) {
			k = rand() % (j+1);
			uiTmp = uiPhi[i][j];
			uiPhi[i][j] = uiPhi[i][k];
			uiPhi[i][k] = uiTmp;
		}
	}

	#ifdef _DEBUG_CFTKD
	DebugOut("Phi permutation (each column)", idx);
	for (int i=0 ; i<iNum ; i++) {
		printf("[DH-%03d] Phi permutation (each column)\t i = %d : ", idx, i);
		for (int j=0 ; j<CLASS_SIZE ; j++)
			printf(" %d ", uiPhi[i][j]);
		std::cout << std::endl;
	}
	#endif

	// ** Computation ** //

	for (int i=0 ; i<iNum ; i++) {
		iTidx = tOut->iRan[idx]+i;

		// tmp1'_i = SM(d'_(i,m+1), TopKRes'_i)
		//SecMul(&cTmp1, &(tIn->cClass[iTidx]), &(tOut->cRes[iTidx]), idx, tRecv, tSend, ucSendPtr);
		SecMul(&cTmp1, &(tIn->cClass[iTidx]), &(tOut->cRes[iTidx]), idx, tRecv);
		#ifdef _DEBUG_CFTKD
		//DebugOut("tmp1'_i = SM(d'_(i,m+1), TopKRes'_i)", cTmp1.c, idx);
		DebugDec("tmp1'_i = SM(d'_(i,m+1), TopKRes'_i)", &cTmp1, idx);
		#endif

		for (int j=0 ; j<CLASS_SIZE ; j++) {
			// tmp2'_i,j = tmp1'_i * (-c_j)' for i=1~n, j=1~v
			paillier_mul(mPubKey, &cT[i][j], &cTmp1, &(tPre->cCls[j]));
			#ifdef _DEBUG_CFTKD
			//DebugOut("tmp2'_i,j = tmp1'_i * (-c_j)'\t", cT[i][j].c, idx);
			DebugDec("tmp2'_i,j = tmp1'_i * (-c_j)'\t", &cT[i][j], idx);
			#endif

			// r'_i,j for i=1~n, j=1~v
			mpz_urandomm(pR.m, state, mPubKey->n);
			#ifdef _DEBUG_CFTKD
			DebugOut("r'_i,j\t\t\t\t", pR.m, idx);
			#endif

			// e'_i,j = tmp2'_i,j ^ r_i,j
			paillier_exp(mPubKey, &cT[i][j], &cT[i][j], &pR);
			#ifdef _DEBUG_CFTKD
			//DebugOut("e'_i,j = tmp2'_i,j ^ r_i,j\t", cT[i][j].c, idx);
			DebugDec("e'_i,j = tmp2'_i,j ^ r_i,j\t", &cT[i][j], idx);
			#endif
		}
	}

	#ifdef _DEBUG_CFTKD
	for (int i=0 ; i<iNum ; i++) {
		for (int j=0 ; j<CLASS_SIZE ; j++) {
			//DebugOut("e'_i,j = tmp2'_i,j ^ r_i,j\t", cT[i][j].c, idx);
			printf("[DH-%03d] (i, j) = (%d, %d) : ", idx, i, j);
			DebugDec(" \t", &cT[i][j], idx);
		}
	}
	#endif

	// thread사이의 데이터 permutation은 살짝 생략하자.
	for (int i=0 ; i<iNum ; i++) {
		k = uiSigma[i];
		for (int j=0 ; j<CLASS_SIZE ; j++) {
			l = uiPhi[k][j];

			// send U'
			paillier_ciphertext_to_bytes(bU, ENC_SIZE, &cT[k][l]);
			SetSendMsg(ucSendPtr, bU, idx, COM_CFTKD, ENC_SIZE);
			#ifdef _DEBUG_CFTKD
		    DebugCom("U' (Hex):\t\t\t", ucSendPtr, ENC_SIZE+HED_SIZE, idx);
			#endif

//		    mSendMtx->lock();
//		    mSendQueue->push(ucSendPtr);
//		    mSendMtx->unlock();
//		    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?

			try {
				mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
			}
			catch ( SocketException& e ) {
				std::cout << "Exception: " << e.description() << std::endl;
			}

			ulRecv.lock();
			tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
			ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
			ulRecv.unlock();

			// W'
			sLen = Byte2Short(ucRecvPtr+HED_LEN);
			#ifdef _DEBUG_CFTKD
		    DebugCom("W' (Hex)\t\t\t", ucRecvPtr, sLen+HED_SIZE, idx);
			#endif

			#ifdef _DEBUG_Assert
		    assert((*(ucRecvPtr+2)==COM_CFTKD)&&(sLen==ENC_SIZE));
			#endif

			paillier_ciphertext_from_bytes(&cW, ucRecvPtr+HED_SIZE, ENC_SIZE);
			#ifdef _DEBUG_CFTKD
			//DebugOut("W (Copied)\t\t\t", cW.c, idx);
			DebugDec("W (Copied)\t\t\t", &cW, idx);
			#endif

			// f'_j = SUM_i=1~n w'_i,j for j=1~v
			paillier_mul(mPubKey, &(tOut->cTF[idx][l]), &(tOut->cTF[idx][l]), &cW);
			#ifdef _DEBUG_CFTKD
			//DebugOut("f'_j = SUM_i=1~n w'_i,j\t", tOut->cTF[idx][l].c, idx);
			DebugDec("f'_j = SUM_i=1~n w'_i,j\t", &(tOut->cTF[idx][l]), idx);
			#endif

		    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
			ucRecvPtr[1] = 0xff;
		}
	}

	#ifdef _DEBUG_CFTKD
	for (int i=0 ; i<CLASS_SIZE ; i++) {
		//DebugOut("i = ", tOut->cTF[i]->c, idx);
		printf("[DH-%03d] i = %d : ", idx, i);
		DebugDec(" \t\t", &(tOut->cTF[idx][i]), idx);
	}
	#endif

	// checking the allocated size  of GMP
	#ifdef _DEBUG_Assert
//	assert(	   pTmp.m->_mp_alloc == 1);
	assert(		 pR.m->_mp_alloc <=   GMP_N_SIZE);
	assert(	  cTmp1.c->_mp_alloc == 2*GMP_N_SIZE*2);
	assert(		 cW.c->_mp_alloc == 2*GMP_N_SIZE);
	assert(cT[0][0].c->_mp_alloc == 2*GMP_N_SIZE*2);
	#endif

	return;
}

void PaillierCrypto::ComputeMc(out_t* tOut, in_t *tIn, unsigned short idx)
{
	paillier_plaintext_t pC;

	#ifdef _DEBUG_INIT_1
	mpz_init(pC.m);
	#else
	mpz_init2(pC.m, 1);
	#endif

	for (int i=tOut->iRan2[idx] ; i<tOut->iRan2[idx+1] ; i++) {
		mpz_set_ui(pC.m, tIn->Class[i]);
		// mc'_j = mf'_j^(c_j)
		paillier_exp(mPubKey, &(tOut->cMc[i]), &(tOut->cFRes[i]), &pC);
	}

	#ifdef _DEBUG_Assert
	assert(pC.m->_mp_alloc == 1);
	#endif

	return;
}

//int PaillierCrypto::TerminatePgm(unsigned short idx, thp_t* tRecv, th_t* tSend, unsigned char* ucSendPtr)
int PaillierCrypto::TerminatePgm(unsigned short idx, thp_t* tRecv)
{
	std::unique_lock<std::mutex> ulRecv(tRecv->m, std::defer_lock);
	unsigned char* ucRecvPtr=NULL;
	unsigned char ucSendPtr[HED_SIZE+1]={0,};
	unsigned char bA[1]={0,};
	short sLen, sRes;

	// send a'
	bA[0] = 0xFF;
	SetSendMsg(ucSendPtr, bA, idx, COM_TERM, 1);
	#ifdef _DEBUG_TERMINATE_PGM
    DebugCom("Terminate Program (sending)\t", ucSendPtr, 1+HED_SIZE, idx);
	#endif

//    mSendMtx->lock();
//    mSendQueue->push(ucSendPtr);
//    mSendMtx->unlock();
//    mSendCV->notify_all();		// Communication thread가 아니라 PPkNN thread가 notify되는 것은 아닌지?

	try {
		mSSocket->send(ucSendPtr, sizeof(ucSendPtr));
	}
	catch ( SocketException& e ) {
		std::cout << "Exception: " << e.description() << std::endl;
	}

	ulRecv.lock();
	tRecv->cv.wait(ulRecv, [&] {return tRecv->pa[idx] != NULL;});
	ucRecvPtr = tRecv->pa[idx];	tRecv->pa[idx] = NULL;
	ulRecv.unlock();

	// received message
	sLen = Byte2Short(ucRecvPtr+HED_LEN);
	#ifdef _DEBUG_TERMINATE_PGM
    DebugCom("Terminate Program (received)\t", ucRecvPtr, sLen+HED_SIZE, idx);
	#endif

	// a' = E(a+r_a)
	#ifdef _DEBUG_Assert
    assert((*(ucRecvPtr+2)==COM_TERM)&&(sLen==3));
	#endif

    sRes = ucRecvPtr[HED_SIZE+2];

    memset(ucRecvPtr, 0, HED_SIZE+ENC_SIZE);
	ucRecvPtr[1] = 0xff;
	return sRes;
}


inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data,
		const unsigned short Idx, const unsigned char Tag, const unsigned short Len)
{
	//memset(ucSendPtr, 0, HED_SIZE+2*ENC_SIZE);
	ucSendPtr[0] = (unsigned char)(Idx & 0x000000ff);
	ucSendPtr[1] = (unsigned char)((Idx & 0x0000ff00) >> 8);
	ucSendPtr[2] = Tag;
	ucSendPtr[3] = (unsigned char)(Len & 0x000000ff);
	ucSendPtr[4] = (unsigned char)((Len & 0x0000ff00) >> 8);
	memcpy(ucSendPtr+HED_SIZE, Data, Len);
	//delete[] Data;

	return;
}

inline void PaillierCrypto::SetSendMsg(unsigned char* ucSendPtr, unsigned char* Data1, unsigned char* Data2,
		const unsigned short Idx, const unsigned char Tag, const unsigned short Len1, const unsigned short Len2)
{
	//memset(ucSendPtr, 0, HED_SIZE+2*ENC_SIZE);
	ucSendPtr[0] = (unsigned char)(Idx & 0x000000ff);
	ucSendPtr[1] = (unsigned char)((Idx & 0x0000ff00) >> 8);
	ucSendPtr[2] = Tag;
	ucSendPtr[3] = (unsigned char)((Len1+Len2) & 0x000000ff);
	ucSendPtr[4] = (unsigned char)(((Len1+Len2) & 0x0000ff00) >> 8);
	memcpy(ucSendPtr+HED_SIZE, Data1, Len1);
	memcpy(ucSendPtr+HED_SIZE+Len1, Data2, Len2);
	//delete[] Data1;
	//delete[] Data2;

	return;
}

//inline void PaillierCrypto::Short2Byte(unsigned char* Ret, unsigned short In)
//{
//	Ret[1] = (unsigned char)((In & 0x0000ff00) >> 8);
//	Ret[0] = (unsigned char)(In & 0x000000ff);
//	return;
//}

inline unsigned short PaillierCrypto::Byte2Short(unsigned char* In)
{
	return (In[1]<<8) + In[0];
}

inline void PaillierCrypto::DebugOut(const mpz_t pData)
{
	gmp_printf("***DD*** %ZX\n", pData);
	return;
}

inline void PaillierCrypto::DebugOut(const char* pMsg, const int idx)
{
	printf("[DH-%03d] %s \n", idx, pMsg);
	return;
}

void PaillierCrypto::DebugOutMain(const char* pMsg, const int idx)
{
	printf("[DH-%03d] %s \n", idx, pMsg);
	return;
}

inline void PaillierCrypto::DebugOut(const char* pMsg, const mpz_t pData)
{
	gmp_printf("2-[DH] %s : %ZX\n", pMsg, pData);
	return;
}

inline void PaillierCrypto::DebugOut(const char* pMsg, const mpz_t pData, const int idx)
{
	gmp_printf("[DH-%03d] %s : %ZX\n", idx, pMsg, pData);
	return;
}

void PaillierCrypto::DebugOutMain(const char* pMsg, const mpz_t pData, const int idx)
{
	gmp_printf("[DH-%03d] %s : %ZX\n", idx, pMsg, pData);
	return;
}

inline void PaillierCrypto::DebugOut(const char* pMsg, const unsigned char* pData)
{
	printf("1-[DH] %s : %s \n", pMsg, pData);
	return;
}

void PaillierCrypto::DebugOut(const char* pMsg, short pData, const short idx)
{
	printf("[DH-%03d] %s : %d\n", idx, pMsg, pData);
	return;
}

inline void PaillierCrypto::DebugOut(paillier_pubkey_t* pPubKey)
{
	gmp_printf("***DD*** (bits, n) %d, %ZX\n", pPubKey->bits, pPubKey->n);
	gmp_printf("***DD*** (n_plusone) %ZX\n", pPubKey->n_plusone);
	gmp_printf("***DD*** (n_squared) %ZX\n", pPubKey->n_squared);
	return;
}

inline void PaillierCrypto::DebugDec(paillier_ciphertext_t* cTxt)
{
	paillier_plaintext_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 3*GMP_N_SIZE);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("5-*DH* %ZX\n", &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc == 3*GMP_N_SIZE);
	#endif

    return;
}

inline void PaillierCrypto::DebugDec(const char* pMsg, paillier_ciphertext_t* cTxt)
{
	paillier_plaintext_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 3*GMP_N_SIZE);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("[DH] %s : %ZX\n", pMsg, &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc == 3*GMP_N_SIZE);
	#endif

    return;
}

inline void PaillierCrypto::DebugDec(const char* pMsg, paillier_ciphertext_t* cTxt, const int idx)
{
	paillier_plaintext_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("[DH-%03d] %s : %ZX\n", idx, pMsg, &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <= 2*GMP_N_SIZE+1);
	#endif

    return;
}

void PaillierCrypto::DebugDecMain(const char* pMsg, paillier_ciphertext_t* cTxt, const int idx)
{
	paillier_plaintext_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("[DH-%03d] %s : %ZX\n", idx, pMsg, &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <=2*GMP_N_SIZE+1);
	#endif

    return;
}

inline void PaillierCrypto::DebugDecBit(paillier_ciphertext_t* cTxt)
{
	paillier_plaintext_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 2*GMP_N_SIZE+1);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("%ZX", &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <= 2*GMP_N_SIZE+1);
	#endif

    return;
}

void PaillierCrypto::DebugDecBitMain(paillier_ciphertext_t* cTxt)
{
	paillier_plaintext_t m;
	#ifdef _DEBUG_INIT_1
	mpz_init(m.m);
	#else
	mpz_init2(m.m, 3*GMP_N_SIZE);
	#endif
	paillier_dec(&m, mPubKey, mPrvKey, cTxt);
	gmp_printf("%ZX", &m);

	#ifdef _DEBUG_Assert
	assert(m.m->_mp_alloc <= 3*GMP_N_SIZE);
	#endif

    return;
}

inline void PaillierCrypto::DebugCom(const char* pMsg, char *data, int len, const int idx)
{
	constexpr char hexmap[] = {'0', '1', '2', '3', '4', '5', '6', '7',
							   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
	char s[len*2+1]={0,};
	for (int i = 0; i < len; ++i) {
		s[2*i]	 = hexmap[(data[i] & 0xF0) >> 4];
		s[2*i+1] = hexmap[data[i] & 0x0F];
	}
	printf("[DH-%03d] %s : %s\n", idx, pMsg, s);

	return;
}

inline void PaillierCrypto::DebugCom(const char* pMsg, unsigned char *data, int len, const int idx)
{
	constexpr char hexmap[] = {'0', '1', '2', '3', '4', '5', '6', '7',
							   '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
	char s[len*2+1]={0,};
	for (int i = 0; i < len; ++i) {
		s[2*i]	 = hexmap[(data[i] & 0xF0) >> 4];
		s[2*i+1] = hexmap[data[i] & 0x0F];
	}
	printf("[DH-%03d] %s : %s\n", idx, pMsg, s);

	return;
}

inline void PaillierCrypto::paillier_ciphertext_from_bytes(paillier_ciphertext_t* ct, void* c, int len )
{
	mpz_import(ct->c, len, 1, 1, 0, 0, c);
	return;
}

inline void PaillierCrypto::paillier_ciphertext_to_bytes(unsigned char* buf, int len, paillier_ciphertext_t* ct )
{
	int cur_len;
	cur_len = mpz_sizeinbase(ct->c, 2);
	cur_len = PAILLIER_BITS_TO_BYTES(cur_len);
	mpz_export(buf + (len - cur_len), 0, 1, 1, 0, 0, ct->c);
	return;
}

//int PaillierCrypto::Sender(th_t* tSend)
//{
//	unsigned char* ucSendPtr;
//	short sIdx, sLen;
//	std::unique_lock<std::mutex> ulSend(tSend->m, std::defer_lock);
//
//	while (1) {
//		ulSend.lock();
//		tSend->cv.wait(ulSend, [&] {return !tSend->q.empty();});
//
//		ucSendPtr = tSend->q.front();
//		tSend->q.pop();
//		ulSend.unlock();
//
//		sIdx = Byte2Short(ucSendPtr);
//		sLen = Byte2Short(ucSendPtr+HED_LEN)+HED_SIZE;
//
//		#ifdef _DEBUG_Communication
//	    DebugCom("(Sender Thread) SendBuf (Hex)\t", ucSendPtr, sLen, sIdx);
//		#endif
//
//		try {
//			mSSocket->send(ucSendPtr, sLen);
//		}
//		catch ( SocketException& e ) {
//			std::cout << "Exception: " << e.description() << std::endl;
//		}
//
//		if (*(ucSendPtr+2) == COM_TERM) {
//			assert(*(ucSendPtr+HED_SIZE) == 0xFF);
//			return 0;
//		}
//
//		memset(ucSendPtr, 0, HED_SIZE+2*ENC_SIZE);
//	}
//}

int PaillierCrypto::Receiver(thp_t* tRecv)
{
	unsigned char* ucRecvBuf[NUM_MAIN_THREAD+1];
	short sIdx, sLen, sbx=0;
	int iReceivedLen;

	for (int i=0 ; i<NUM_MAIN_THREAD+1 ; i++) {
		ucRecvBuf[i] = new unsigned char[HED_SIZE+ENC_SIZE];
		memset(ucRecvBuf[i], 0, HED_SIZE+ENC_SIZE);
		ucRecvBuf[i][1] = 0xff;
	}

	while (1) {
		while (ucRecvBuf[sbx][1]!=0xff) {
			sbx++;
			if (sbx>=NUM_MAIN_THREAD+1)
				sbx=0;
		}
		iReceivedLen = 0;

		//?? exception 출력, connection reset by peer (errno=104) 처리
		try {
			//while (memcmp(ucRecvBuf[sbx]+1, ucTemp+1, HED_SIZE-1)==0) {
			//	iReceivedLen = mSSocket->recv(ucRecvBuf[sbx], HED_SIZE);
			//}
			iReceivedLen = mSSocket->recv(ucRecvBuf[sbx], HED_SIZE);
			sLen = Byte2Short(ucRecvBuf[sbx]+HED_LEN);
			while (iReceivedLen != HED_SIZE+sLen)
				iReceivedLen += mSSocket->recv(ucRecvBuf[sbx]+iReceivedLen, HED_SIZE+sLen-iReceivedLen);
	    }
		catch ( SocketException& e ) {
			std::cout << "Exception: " << e.description() << std::endl;
		}

		sIdx = Byte2Short(ucRecvBuf[sbx]);

		#ifdef _DEBUG_Communication
		DebugCom("(Receiver Thread) RecvBuf (Hex)", ucRecvBuf[sbx], HED_SIZE+sLen, sIdx);
		#endif

		if (*(ucRecvBuf[sbx]+2) == COM_TERM) {
			assert(Byte2Short(ucRecvBuf[sbx]+HED_SIZE) == MOD_SIZE);
			unsigned char* ucKillMsg = new unsigned char[HED_SIZE+3];
			memcpy(ucKillMsg, ucRecvBuf[sbx], HED_SIZE+3);
		    tRecv->m.lock();
		    tRecv->pa[sIdx] = ucKillMsg;
		    tRecv->cv.notify_all();
		    tRecv->m.unlock();

			for (int i=0 ; i<NUM_MAIN_THREAD ; i++)
				delete[] ucRecvBuf[i];
			return 0;
		}

		tRecv->m.lock();
	    tRecv->pa[sIdx] = ucRecvBuf[sbx];
	    tRecv->m.unlock();
	    tRecv->cv.notify_all();
	}

	return 0;
}
