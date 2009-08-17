#include <Python.h>
#include <string>
#include <gmp.h>
#include <math.h>
#include <ctime> 
#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>


static PyObject *rand(PyObject *self, PyObject *args){
    unsigned int seed = 123456;

    mpz_t rand;
    mpz_init(rand);
    gmp_randstate_t rand_state;

    gmp_randinit_mt(rand_state);
    gmp_randseed_ui(rand_state, seed);

    mpz_urandomb(rand, rand_state, 16);
    char *randout = new char[mpz_sizeinbase(rand, 10)+2];
    mpz_get_str(randout, 10, rand);

    gmp_randclear(rand_state);
    mpz_clear(rand);
    return Py_BuildValue("s", randout);
}
    
static PyObject *generateKeys(PyObject *self, PyObject *args){
    unsigned long int bits;
    if (!PyArg_ParseTuple(args, "K", &bits))
        return NULL;
    unsigned int keyLength = (bits + 3)/4;


    ///////////////////////////////////////
    // initialize!
    gmp_randstate_t rand_state;    
    mpz_t largePrime_P;
    mpz_t largePrime_Q;
    mpz_t coPrime_E;
    mpz_t GCD;
    mpz_t seed;
    mpz_t random;
    mpz_t rand_add;
    mpz_t modulus;
    mpz_t length;
    mpz_t keyN;
    mpz_t keyD;
    mpz_t t1;
    mpz_t t2;
    mpz_t t3;
    mpz_t keyM;
    
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(t3);
    mpz_init(keyM);
    mpz_init(keyN);
    mpz_init(keyD);
    mpz_init(coPrime_E);
    mpz_init(modulus);
    mpz_init(GCD);
    mpz_init(random);
    mpz_init(length);
    mpz_init(rand_add);

    // Setup the random state
    time_t seconds;
    seconds = time (NULL);

    ///////////////////////////////////////
    // setup the seed for random
    gmp_randinit_default(rand_state);
    mpz_init_set_ui(seed, seconds);
    mpz_mul_ui(seed, seed, 31415926);
    mpz_mul_ui(seed, seed, keyLength);
    gmp_randseed(rand_state, seed);

    ///////////////////////////////////////
    // randomly get large prime P
    mpz_init_set_ui(largePrime_P, 0);
    int p = keyLength;
	while ((p -= 2) > 0){
        mpz_mul_ui(largePrime_P, largePrime_P, 16);
        mpz_urandomb(random, rand_state, 8);
		mpz_add(largePrime_P, largePrime_P, random);
	}

    ///////////////////////////////////////
    // make sure its a prime
    while (!mpz_probab_prime_p(largePrime_P, 25)){
        mpz_add_ui(largePrime_P, largePrime_P, 1);
    }

    ///////////////////////////////////////
    // randomly get large prime Q
    mpz_init_set_ui(largePrime_Q, 0);
    int q = keyLength;
    while ((q -= 2) > 0){
        mpz_mul_ui(largePrime_Q, largePrime_Q, 16);
        mpz_urandomb(rand_add, rand_state, 8);
        // urandomb(..., 8) 2^n-1 = 255
        mpz_add(largePrime_Q, largePrime_Q, rand_add);

    }

    ///////////////////////////////////////
    // make sure its a prime
    while (!mpz_probab_prime_p(largePrime_Q, 25)){
        mpz_add_ui(largePrime_Q, largePrime_Q, 1);
    }

    ///////////////////////////////////////
    // multiply P and Q into N
    mpz_mul(keyN, largePrime_P, largePrime_Q);
    char *key_N = new char[mpz_sizeinbase(keyN, 10)+2];
    mpz_get_str(key_N, 10, keyN);
    
    ///////////////////////////////////////
    // M = (P-1) * (Q-1) 
    mpz_sub_ui(t1, largePrime_P, 1);
    mpz_sub_ui(t2, largePrime_Q, 1);
    mpz_mul(keyM, t1, t2);

    ///////////////////////////////////////
    // Coprime E to M
    mpz_init_set_ui(coPrime_E, 0);
    int z = keyLength;
    while ((z--) > 0){
        mpz_mul_ui(coPrime_E, coPrime_E, 16);
        mpz_urandomb(rand_add, rand_state, 8);
        // urandomb(..., 8) 2^n-1 = 255
        mpz_add(coPrime_E, coPrime_E, rand_add);

    }
    while (mpz_cmp(coPrime_E, keyN) >= 0)
			mpz_div_ui(coPrime_E, coPrime_E, 2);

	while (mpz_cmp_ui(GCD, 1) != 1 && mpz_cmp(GCD, keyM) <= 0){ 
	    mpz_sub_ui(coPrime_E, coPrime_E, 1);
        mpz_gcd(GCD, coPrime_E, keyM);
    }
    char *coPrimeE = new char[mpz_sizeinbase(coPrime_E, 10)+2];
    mpz_get_str(coPrimeE, 10, coPrime_E);
    
    ///////////////////////////////////////
    // Find D such that DE % M = 1
    mpz_init_set_ui(keyD, 0);		
    z = keyLength;
    while ((z--) > 0){
        mpz_mul_ui(keyD, keyD, 16);
        mpz_urandomb(rand_add, rand_state, 8);
        // urandomb(..., 8) 2^n-1 = 255
        mpz_add(keyD, keyD, rand_add);
    }
    while (mpz_cmp_ui(modulus,1) != 1){
        mpz_sub_ui(keyD, keyD, 1);
        mpz_mul(t3, keyD, coPrime_E);
        mpz_mod(modulus, t3, keyM);
    }

    char *key_D = new char[mpz_sizeinbase(keyD, 10)+2];
    mpz_get_str(key_D, 10, keyD);


    mpz_clear(t1);
    mpz_clear(t2);
    mpz_clear(t3);
    mpz_clear(keyM);
    mpz_clear(keyN);
    mpz_clear(keyD);
    mpz_clear(coPrime_E);
    mpz_clear(modulus);
    mpz_clear(GCD);
    mpz_clear(random);
    mpz_clear(length);
    mpz_clear(rand_add);

    // E = exponent
    // N = modulus
    // D = exponent
    // N = modulus
    return Py_BuildValue("(ssss)", key_N, coPrimeE, key_N, key_D);
}


static PyObject *encrypt(PyObject *self, PyObject *args){
    const char* toBeEnc;
    const char* publicKey_modulus;
    const char* publicKey_exponent;
    if (!PyArg_ParseTuple(args, "s(ss)", &toBeEnc, &publicKey_modulus, &publicKey_exponent))
        return NULL;


    ///////////////////////////////////////
    // Initialization
    std::string modulus = publicKey_modulus;
    std::string encrypting  = toBeEnc; 
    const unsigned long int chunkSize(((modulus.length() - 2) / 3));
    const unsigned long int chunkCount = encrypting.length() / chunkSize;

    mpz_t encrypted;
    mpz_t toEncrypt;
    mpz_t keyE;
    mpz_t keyN;
    mpz_init(keyE);
    mpz_init(keyN);
    mpz_init(toEncrypt);
    mpz_init(encrypted);
   
    mpz_set_str(keyN, publicKey_modulus, 10);
    mpz_set_str(keyE, publicKey_exponent, 10);

    std::string encoded;
    std::string cypherText;

    ///////////////////////////////////////
    // Breakup into chunks
    for (unsigned long int i(0); i < chunkCount; i++){
        std::string chunk(encrypting.substr(i * chunkSize, chunkSize));

        ///////////////////////////////////////
        // Convert to decimal
        encoded.clear();
        encoded.resize(encrypting.length() * 3 + 1);
        unsigned long int index = encrypting.length() * 3;
        // Factor this out to two methods encode/decode a2i/i2a
        for (unsigned long int k(0); k < encrypting.length(); k++){
                // Encode the characters using their ASCII values' digits
                unsigned char ASCII = encrypting[k];
                encoded[index - 2] = (ASCII % 10) + '0';
                ASCII /= 10;
                encoded[index - 1] = (ASCII % 10) + '0';
                encoded[index] = (ASCII / 10) + '0';
                index -= 3;
        }
        // Convert into a big num for gmp
        mpz_set_str(toEncrypt, encoded.c_str(), 10);
        encoded[0] = '1';

        ///////////////////////////////////////
        // Encrypt
        mpz_powm(encrypted, toEncrypt, keyE, keyN);
        char *encryptedChunk = new char[mpz_sizeinbase(encrypted, 10)+2];
        mpz_get_str(encryptedChunk, 10, encrypted);
        cypherText.append(encryptedChunk);
        cypherText.append(" ");
        delete encryptedChunk;
    }
    ///////////////////////////////////////
    // the message is shorter than the chunk
    if(chunkCount == 0){
        encoded.resize(encrypting.length() * 3 + 1);
        unsigned long int index = encrypting.length() * 3;
        for (unsigned long int k(0); k < encrypting.length(); k++){
                // Encode the characters using their ASCII values' digits
                unsigned char ASCII = encrypting[k];
                encoded[index - 2] = (ASCII % 10) + '0';
                ASCII /= 10;
                encoded[index - 1] = (ASCII % 10) + '0';
                encoded[index] = (ASCII / 10) + '0';
                index -= 3;
        }
        encoded[0] = '1';
        mpz_set_str(toEncrypt, encoded.c_str(), 10);
        ///////////////////////////////////////
        // Encrypt
        mpz_powm(encrypted, toEncrypt, keyE, keyN);
        char *encryptedChunk = new char[mpz_sizeinbase(encrypted, 10)+2];
        mpz_get_str(encryptedChunk, 10, encrypted);
        cypherText = encryptedChunk;
        cypherText.append(" ");
        delete encryptedChunk;
    }

 
    ///////////////////////////////////////
    // Cleanup
    mpz_clear(keyE);
    mpz_clear(keyN);
    mpz_clear(encrypted);
    mpz_clear(toEncrypt);

    return Py_BuildValue("s", cypherText.c_str());
}

static PyObject *decrypt(PyObject *self, PyObject *args){
    const char* toBeDec;
    const char* privateKey_modulus;
    const char* privateKey_exponent;
    if (!PyArg_ParseTuple(args, "s(ss)", &toBeDec, &privateKey_modulus, &privateKey_exponent))
        return NULL;

    ///////////////////////////////////////
    // Initialize
    std::string modulus     = privateKey_modulus;
    std::string cypherText  = toBeDec; 
    std::string message;
    const unsigned long int chunkSize(((modulus.length() - 2) / 3));
    const unsigned long int chunkCount = cypherText.length() / chunkSize;

    mpz_t decrypted;
    mpz_t toDecrypt;
    mpz_t keyD;
    mpz_t keyN;
    mpz_init(keyD);
    mpz_init(keyN);
    mpz_init(toDecrypt);
    mpz_init(decrypted);

    mpz_set_str(keyN, privateKey_modulus, 10);
    mpz_set_str(keyD, privateKey_exponent, 10);
    mpz_set_str(toDecrypt, toBeDec, 10);

    ///////////////////////////////////////
    // Dechunk
    long int i(0), j(0);
    if(chunkCount >= 1){
        while ((j = cypherText.find(' ', i)) != -1){
            std::string chunk = cypherText.substr(i, j - i);
            if (chunk.length() >= modulus.length()){
                PyErr_SetString(PyExc_RuntimeError, "chunk too big");
            }
            mpz_set_str(toDecrypt, chunk.c_str(), 10);
            mpz_powm(decrypted, toDecrypt, keyD, keyN);
            char *decryptedChunk = new char[mpz_sizeinbase(decrypted, 10)+2];
            mpz_get_str(decryptedChunk, 10, decrypted);
            message.append(decryptedChunk);
            delete decryptedChunk;
            i = j + 1;

        }
    } else{
        mpz_powm(decrypted, toDecrypt, keyD, keyN);
        char *decryptedChunk = new char[mpz_sizeinbase(decrypted, 10)+2];
        mpz_get_str(decryptedChunk, 10, decrypted);
        message.append(decryptedChunk);
        delete decryptedChunk;
    }

    
    ///////////////////////////////////////
    // Decode from numerical to ASCII
    std::string decoded;
    for (unsigned long int i(0); i < message.length() / 3; i++){
        // Decode the characters using the ASCII values
        char ASCII = 100 * char(message[i * 3]);
        ASCII += 10 * char(message[i * 3 + 1]);
        decoded.push_back(ASCII + char(message[i * 3 + 2]));
    }
    
    mpz_clear(keyD);
    mpz_clear(keyN);
    mpz_clear(decrypted);
    mpz_clear(toDecrypt);
    return Py_BuildValue("ss", decoded.c_str(), message.c_str());
}


static PyMethodDef RSA_Methods[] = {
    {"encrypt", (PyCFunction)encrypt, METH_VARARGS, ""},
    {"decrypt", (PyCFunction)decrypt, METH_VARARGS, ""},
    {"genKeys", (PyCFunction)generateKeys, METH_VARARGS, ""},
    {"rand", (PyCFunction)rand, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

PyMODINIT_FUNC initrsa(void){
  /* Create the module and add the functions */
  (void)Py_InitModule("rsa", RSA_Methods);
}
