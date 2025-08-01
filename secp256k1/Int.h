#ifndef BIGINTH
#define BIGINTH

#include <string>
#include <inttypes.h>
#include <gmp.h>

// We need 1 extra block for Knuth div algorithm , Montgomery multiplication and ModInv
#define BISIZE 256

#if BISIZE==256
  #define NB64BLOCK 5
  #define NB32BLOCK 10
#elif BISIZE==512
  #define NB64BLOCK 9
  #define NB32BLOCK 18
#else
  #error Unsuported size
#endif

class Int {

public:

  Int();
  Int(int64_t i64);
  Int(uint64_t u64);
  Int(Int *a);

  // Op
  void floor_Div(Int *a, Int *b);
  void Add(uint64_t a);
  void Add(Int *a);
  void Add(Int *a, Int *b);
  void AddOne();
  void Sub(uint64_t a);
  void Sub(Int *a);
  void Sub(Int *a, Int *b);
  void SubOne();
  void Mult(Int *a);
  void Mult(uint64_t a);
  void IMult(int64_t a);
  void Mult(Int *a, uint64_t b);
  void IMult(Int *a, int64_t b);
  void Mult(Int *a, Int *b);
  void IMultAdd(Int* a, int64_t aa, Int* b, int64_t bb);
  void Div(Int *a, Int *mod = NULL);
  void MultModN(Int *a, Int *b, Int *n);
  void Neg();
  void Abs();

  // Right shift (signed)
  void ShiftR(uint32_t n);
  void ShiftR32Bit();
  void ShiftR64Bit();
  // Left shift
  void ShiftL(uint32_t n);
  void ShiftL32Bit();
  void ShiftL64Bit();
  // Bit swap
  void SwapBit(int bitNumber);


  // Comp
  bool IsGreater(Int *a);
  bool IsGreaterOrEqual(Int *a);
  bool IsLowerOrEqual(Int *a);
  bool IsLower(Int *a);
  bool IsEqual(Int *a);
  bool IsZero();
  bool IsOne();
  bool IsStrictPositive();
  bool IsPositive();
  bool IsNegative();
  bool IsEven();
  bool IsOdd();

  double ToDouble();

  // Modular arithmetic

  // Setup field
  // n is the field characteristic
  // R used in Montgomery mult (R = 2^size(n))
  // R2 = R^2, R3 = R^3, R4 = R^4
  static void SetupField(Int *n, Int *R = NULL, Int *R2 = NULL, Int *R3 = NULL, Int *R4 = NULL);
  static Int *GetR();                            // Return R
  static Int *GetR2();                           // Return R2
  static Int *GetR3();                           // Return R3
  static Int *GetR4();                           // Return R4
  static Int* GetFieldCharacteristic();          // Return field characteristic

  void GCD(Int *a);                          // this <- GCD(this,a)
  void Mod(Int *n);                          // this <- this (mod n)
  void ModInv();                             // this <- this^-1 (mod n)
  void MontgomeryMult(Int *a,Int *b);        // this <- a*b*R^-1 (mod n)
  void MontgomeryMult(Int *a);               // this <- this*a*R^-1 (mod n)
  void ModAdd(Int *a);                       // this <- this+a (mod n) [0<a<P]
  void ModAdd(Int *a,Int *b);                // this <- a+b (mod n) [0<a,b<P]
  void ModAdd(uint64_t a);                   // this <- this+a (mod n) [0<a<P]
  void ModSub(Int *a);                       // this <- this-a (mod n) [0<a<P]
  void ModSub(Int *a, Int *b);               // this <- a-b (mod n) [0<a,b<P]
  void ModSub(uint64_t a);                   // this <- this-a (mod n) [0<a<P]
  void ModMul(Int *a,Int *b);                // this <- a*b (mod n)
  void ModMul(Int *a);                       // this <- this*b (mod n)
  void ModSquare(Int *a);                    // this <- a^2 (mod n)
  void ModCube(Int *a);                      // this <- a^3 (mod n)
  void ModDouble();                          // this <- 2*this (mod n)
  void ModExp(Int *e);                       // this <- this^e (mod n)
  void ModNeg();                             // this <- -this (mod n)
  void ModSqrt();                            // this <- +/-sqrt(this) (mod n)
  bool HasSqrt();                            // true if this admit a square root

  // Specific SecpK1
  static void InitK1(Int *order);
  void ModMulK1(Int *a, Int *b);
  void ModMulK1(Int *a);
  void ModSquareK1(Int *a);
  void ModMulK1order(Int *a);
  void ModAddK1order(Int *a,Int *b);
  void ModAddK1order(Int *a);
  void ModSubK1order(Int *a);
  void ModNegK1order();
  uint32_t ModPositiveK1();

  // Size
  int GetSize();
  int GetBitLength();

  // Setter
  void SetInt32(uint32_t value);
  void SetInt64(uint64_t value);
  void Set(Int *a);
  void SetBase10(char *value);
  void SetBase16(char *value);
  void SetBaseN(int n,char *charset,char *value);
  void SetByte(int n,unsigned char byte);
  void SetDWord(int n, uint32_t b);
  void SetQWord(int n,uint64_t b);
  void Set32Bytes(unsigned char *bytes);
  void MaskByte(int n);

  // Getter
  uint32_t GetInt32();
  int GetBit(uint32_t n);
  unsigned char GetByte(int n);
  void Get32Bytes(unsigned char *buff);

  // To String
  std::string GetBase2();
  std::string GetBase10();
  std::string GetBase16();
  std::string GetBaseN(int n, char *charset);
  std::string GetBlockStr();
  std::string GetC64Str(int nbDigit);

  // Check function
  static void Check();


  /*
  // Align to 16 bytes boundary
  union {
    __declspec(align(16)) uint32_t bits[NB32BLOCK];
    __declspec(align(16)) uint64_t bits64[NB64BLOCK];
  };
  */
  union {
    uint32_t bits[NB32BLOCK];
    uint64_t bits64[NB64BLOCK];
  };

private:

  void ShiftL32BitAndSub(Int *a, int n);
  uint64_t AddC(Int *a);
  void AddAndShift(Int *a, Int *b, uint64_t cH);
  void Mult(Int *a, uint32_t b);
  int  GetLowestBit();
  void CLEAR();
  void CLEARFF();

};

// Inline routines

#ifndef WIN64

// Missing intrinsics
static uint64_t inline _umul128(uint64_t a, uint64_t b, uint64_t *h) {
  uint64_t rhi;
  uint64_t rlo;
  __asm__( "mulq  %[b];" :"=d"(rhi),"=a"(rlo) :"1"(a),[b]"rm"(b));
    *h = rhi;
    return rlo;
}

static uint64_t inline __shiftright128(uint64_t a, uint64_t b, unsigned char n) {
  uint64_t c;
  __asm__ ("movq %1,%0;shrdq %3,%2,%0;" : "=D"(c) : "r"(a),"r"(b),"c"(n));
  return  c;
}


static uint64_t inline __shiftleft128(uint64_t a, uint64_t b, unsigned char n) {
  uint64_t c;
  __asm__ ("movq %1,%0;shldq %3,%2,%0;" : "=D"(c) : "r"(b),"r"(a),"c"(n));
  return  c;
}

#define _subborrow_u64(a,b,c,d) __builtin_ia32_sbb_u64(a,b,c,(long long unsigned int*)d);
//#define _subborrow_u64(a,b,c,d) __builtin_ia32_subborrow_u64(a,b,c,(long long unsigned int*)d); // clang++
#define _addcarry_u64(a,b,c,d) __builtin_ia32_addcarryx_u64(a,b,c,(long long unsigned int*)d);
#define _byteswap_uint64 __builtin_bswap64

#else

#include <intrin.h>

static inline int __builtin_ctzll(unsigned long long x) {
  unsigned long ret;
  _BitScanForward64(&ret,x);
  return (int)ret;
}

#endif


#define LoadI64(i,i64)    \
i.bits64[0] = i64;        \
i.bits64[1] = i64 >> 63;  \
i.bits64[2] = i.bits64[1];\
i.bits64[3] = i.bits64[1];\
i.bits64[4] = i.bits64[1];

static void inline imm_mul(uint64_t *x, uint64_t y, uint64_t *dst) {

  unsigned char c = 0;
  uint64_t h, carry;
  dst[0] = _umul128(x[0], y, &h); carry = h;
  c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
  c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
  c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;
  c = _addcarry_u64(c, _umul128(x[4], y, &h), carry, dst + 4); carry = h;
#if NB64BLOCK > 5
  c = _addcarry_u64(c, _umul128(x[5], y, &h), carry, dst + 5); carry = h;
  c = _addcarry_u64(c, _umul128(x[6], y, &h), carry, dst + 6); carry = h;
  c = _addcarry_u64(c, _umul128(x[7], y, &h), carry, dst + 7); carry = h;
  c = _addcarry_u64(c, _umul128(x[8], y, &h), carry, dst + 8); carry = h;
#endif

}

static void inline imm_umul(uint64_t *x, uint64_t y, uint64_t *dst) {

  // Assume that x[NB64BLOCK-1] is 0
  unsigned char c = 0;
  uint64_t h, carry;
  dst[0] = _umul128(x[0], y, &h); carry = h;
  c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
  c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
  c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;
#if NB64BLOCK > 5
  c = _addcarry_u64(c, _umul128(x[4], y, &h), carry, dst + 4); carry = h;
  c = _addcarry_u64(c, _umul128(x[5], y, &h), carry, dst + 5); carry = h;
  c = _addcarry_u64(c, _umul128(x[6], y, &h), carry, dst + 6); carry = h;
  c = _addcarry_u64(c, _umul128(x[7], y, &h), carry, dst + 7); carry = h;
#endif
  _addcarry_u64(c, 0ULL, carry, dst + (NB64BLOCK - 1));

}

static void inline shiftR(unsigned char n, uint64_t *d) {

  d[0] = __shiftright128(d[0], d[1], n);
  d[1] = __shiftright128(d[1], d[2], n);
  d[2] = __shiftright128(d[2], d[3], n);
  d[3] = __shiftright128(d[3], d[4], n);
#if NB64BLOCK > 5
  d[4] = __shiftright128(d[4], d[5], n);
  d[5] = __shiftright128(d[5], d[6], n);
  d[6] = __shiftright128(d[6], d[7], n);
  d[7] = __shiftright128(d[7], d[8], n);
#endif
  d[NB64BLOCK-1] = ((int64_t)d[NB64BLOCK-1]) >> n;

}

static void inline shiftL(unsigned char n, uint64_t *d) {

#if NB64BLOCK > 5
  d[8] = __shiftleft128(d[7], d[8], n);
  d[7] = __shiftleft128(d[6], d[7], n);
  d[6] = __shiftleft128(d[5], d[6], n);
  d[5] = __shiftleft128(d[4], d[5], n);
#endif
  d[4] = __shiftleft128(d[3], d[4], n);
  d[3] = __shiftleft128(d[2], d[3], n);
  d[2] = __shiftleft128(d[1], d[2], n);
  d[1] = __shiftleft128(d[0], d[1], n);
  d[0] = d[0] << n;

}

#endif // BIGINTH
