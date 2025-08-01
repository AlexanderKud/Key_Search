#include "Int.h"
#include <emmintrin.h>
#include <string.h>

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))

static Int     Field_P;       // Field characteristic
static Int     _R;       // Montgomery multiplication R
static Int     _R2;      // Montgomery multiplication R2
static Int     _R3;      // Montgomery multiplication R3
static Int     _R4;      // Montgomery multiplication R4
static int32_t  Msize;    // Montgomery mult size
static uint32_t MM32;     // 32bits lsb negative inverse of P
static uint64_t MM64;     // 64bits lsb negative inverse of P
#define MSK62  0x3FFFFFFFFFFFFFFF

extern Int _ONE;

// ------------------------------------------------

void Int::ModAdd(Int *a) {
  Int p;
  Add(a);
  p.Sub(this, &Field_P);
  if(p.IsPositive())
    Set(&p);
}

// ------------------------------------------------

void Int::ModAdd(Int *a, Int *b) {
  Int p;
  Add(a,b);
  p.Sub(this, &Field_P);
  if(p.IsPositive())
    Set(&p);
}

// ------------------------------------------------

void Int::ModDouble() {
  Int p;
  Add(this);
  p.Sub(this, &Field_P);
  if(p.IsPositive())
    Set(&p);
}

// ------------------------------------------------

void Int::ModAdd(uint64_t a) {
  Int p;
  Add(a);
  p.Sub(this, &Field_P);
  if(p.IsPositive())
    Set(&p);
}

// ------------------------------------------------

void Int::ModSub(Int *a) {
  Sub(a);
  if (IsNegative())
    Add(&Field_P);
}

// ------------------------------------------------

void Int::ModSub(uint64_t a) {
  Sub(a);
  if (IsNegative())
    Add(&Field_P);
}

// ------------------------------------------------

void Int::ModSub(Int *a, Int *b) {
  Sub(a, b);
  if (IsNegative())
    Add(&Field_P);
}

// ------------------------------------------------

void Int::ModNeg() {
  Neg();
  Add(&Field_P);
}

// ------------------------------------------------

inline void DivStep62(int64_t u0,int64_t v0,
  int64_t* eta,
  int64_t* uu,int64_t* uv,
  int64_t* vu,int64_t* vv) {

  // u' = (uu*u + uv*v) >> bitCount
  // v' = (vu*u + vv*v) >> bitCount
  // Do not maintain a matrix for r and s, the number of
  // 'added P' can be easily calculated

  int bitCount;

#if 0

  #define SWAP_ADD(x,y) x+=y;y-=x;
  #define SWAP_SUB(x,y) x-=y;y+=x;

  // Former divstep62 (using __builtin_ctzll)
  // Do not use eta, u and v have an exponential decay in worst case
  // but with low probability to reach this worst case complexity
  // Avg: 581 Kinv/s

  bitCount = 62;
  int64_t nb0;

  while(true) {

    int zeros = __builtin_ctzll(v0 | (UINT64_MAX << bitCount));
    v0 >>= zeros;
    *uu <<= zeros;
    *uv <<= zeros;
    bitCount -= zeros;

    if(bitCount <= 0)
      break;

    nb0 = (v0 + u0) & 0x3;
    if(nb0 == 0) {
      SWAP_ADD(*vv,*uv);
      SWAP_ADD(*vu,*uu);
      SWAP_ADD(v0,u0);
    } else {
      SWAP_SUB(*vv,*uv);
      SWAP_SUB(*vu,*uu);
      SWAP_SUB(v0,u0);
    }

  }


#endif

#if 1

  #define SWAP_NEG(tmp,x,y) tmp = x; x = y; y = -tmp;

  int64_t m,w,x,y,z;
  bitCount = 62;

  // divstep62 var time implementation by Peter Dettman
  // (see https://github.com/bitcoin-core/secp256k1/pull/767)
  // Avg: 640 Kinv/s

  while(true) {

    // Use a sentinel bit to count zeros only up to bitCount
    int zeros = __builtin_ctzll(v0 | (UINT64_MAX << bitCount));

    v0 >>= zeros;
    *uu <<= zeros;
    *uv <<= zeros;
    *eta -= zeros;
    bitCount -= zeros;

    if(bitCount <= 0)
      break;

    if(*eta < 0) {
      *eta = -*eta;
      SWAP_NEG(x,u0,v0);
      SWAP_NEG(y,*uu,*vu);
      SWAP_NEG(z,*uv,*vv);
    }

    // Handle up to 3 divsteps at once, subject to eta and bitCount
    int limit = (*eta + 1) > bitCount ? bitCount : (int)(*eta + 1);
    m = (UINT64_MAX >> (64 - limit)) & 7U;

    // Note that f * f == 1 mod 8, for any f
    w = (-u0 * v0) & m;
    v0 += u0 * w;
    *vu += *uu * w;
    *vv += *uv * w;

  }

#endif

#if 0

  // divstep62 constant time implementation by Peter Dettman
  // Avg: 405 Kinv/s

  uint64_t c1,c2,x,y,z;

  for(bitCount = 0; bitCount < 62; bitCount++) {

    c1 = -(v0 & ((uint64_t)(*eta) >> 63));

    x = (u0 ^ v0) & c1;
    u0 ^= x; v0 ^= x; v0 ^= c1; v0 -= c1;

    y = (*uu ^ *vu) & c1;
    *uu ^= y; *vu ^= y; *vu ^= c1; *vu -= c1;

    z = (*uv ^ *vv) & c1;
    *uv ^= z; *vv ^= z; *vv ^= c1; *vv -= c1;

    *eta = (*eta ^ c1) - c1 - 1;

    c2 = -(v0 & 1);

    v0 += (u0 & c2); v0 >>= 1;
    *vu += (*uu & c2); *uu <<= 1;
    *vv += (*uv & c2); *uv <<= 1;
  }

#endif

}

#define MatrixVecMul(u,v,_11,_12,_21,_22) \
  t1.IMult(&u,_11); \
  t2.IMult(&v,_12); \
  t3.IMult(&u,_21); \
  t4.IMult(&v,_22); \
  u.Add(&t1,&t2); \
  v.Add(&t3,&t4);

// ------------------------------------------------

void Int::ModInv() {

  // Compute modular inverse of this mop _P
  // 0 < this < P  , P must be odd
  // Return 0 if no inverse

  // 256bit
  //#define XCD 1               // ~80  kOps/s
  //#define MONTGOMERY 1        // ~246 kOps/s
  //#define PENK 1              // ~215 kOps/s
  #define DRS62 1               // ~640 kOps/s

  Int u(&Field_P);
  Int v(this);

#ifdef XCD

  Int r((int64_t)0);
  Int s((int64_t)1);
  Int q, t1, t2, w;

  // Classic XCD

  bool bIterations = true;  // Remember odd/even iterations
  while (!u.IsZero()) {
    // Step X3. Divide and "Subtract"
    q.Set(&v);
    q.Div(&u, &t2);   // q = u / v, t2 = u % v
    w.Mult(&q, &r);   // w = q * r
    t1.Add(&s, &w);   // t1 = s + w

                      // Swap u,v & r,s
    s.Set(&r);
    r.Set(&t1);
    v.Set(&u);
    u.Set(&t2);

    bIterations = !bIterations;
  }

  if (!v.IsOne()) {
    CLEAR();
    return;
  }

  if (!bIterations) {
    Set(&Field_P);
    Sub(&s);  /* inv = n - u1 */
  } else {
    Set(&s);  /* inv = u1     */
  }

#endif

#ifdef PENK

  Int r((int64_t)0);
  Int s((int64_t)1);
  Int x;
  Int n2(&Field_P);
  int k = 0;
  int T;
  int Q = Field_P.bits[0] & 3;
  shiftL(1,n2.bits64);

  // Penk's Algorithm (With DRS2 optimisation)

  while (v.IsEven()) {

    shiftR(1,v.bits64);
    if (s.IsEven())
      shiftR(1, s.bits64);
    else if (s.IsGreater(&Field_P)) {
      s.Sub(&Field_P);
      shiftR(1, s.bits64);
    } else {
      s.Add(&Field_P);
      shiftR(1, s.bits64);
    }

  }

  while (true) {

    if (u.IsGreater(&v)) {

      if ((u.bits[0] & 2) == (v.bits[0] & 2)) {
        u.Sub(&v);
        r.Sub(&s);
      } else {
        u.Add(&v);
        r.Add(&s);
      }
      shiftR(2,u.bits64);
      T = r.bits[0] & 3;
      if (T == 0) {
        shiftR(2,r.bits64);
      } else if (T == 2) {
        r.Add(&n2);
        shiftR(2, r.bits64);
      } else if (T == Q) {
        r.Sub(&Field_P);
        shiftR(2, r.bits64);
      } else {
        r.Add(&Field_P);
        shiftR(2, r.bits64);
      }
      while (u.IsEven()) {
        shiftR(1,u.bits64);
        if (r.IsEven()) {
          shiftR(1, r.bits64);
        } else if (r.IsGreater(&Field_P)) {
          r.Sub(&Field_P);
          shiftR(1, r.bits64);
        } else {
          r.Add(&Field_P);
          shiftR(1, r.bits64);
        }
      }

    } else {

      if ((u.bits[0] & 2) == (v.bits[0] & 2)) {
        v.Sub(&u);
        s.Sub(&r);
      } else {
        v.Add(&u);
        s.Add(&r);
      }

      if (v.IsZero())
        break;

      shiftR(2, v.bits64);
      T = s.bits[0] & 3;
      if (T == 0) {
        shiftR(2,s.bits64);
      } else if (T == 2) {
        s.Add(&n2);
        shiftR(2, s.bits64);
      } else if (T == Q) {
        s.Sub(&Field_P);
        shiftR(2, s.bits64);
      } else {
        s.Add(&Field_P);
        shiftR(2, s.bits64);
      }

      while (v.IsEven()) {
        shiftR(1, v.bits64);
        if (s.IsEven()) {
          shiftR(1, s.bits64);
        } else if (s.IsGreater(&Field_P)) {
          s.Sub(&Field_P);
          shiftR(1, s.bits64);
        } else {
          s.Add(&Field_P);
          shiftR(1, s.bits64);
        }
      }

    }

  }

  if (u.IsGreater(&_ONE)) {
    CLEAR();
    return;
  }
  if (r.IsNegative())
    r.Add(&Field_P);
  Set(&r);

#endif

#ifdef MONTGOMERY

  Int r((int64_t)0);
  Int s((int64_t)1);
  Int x;
  int k = 0;

  // Montgomery method
  while (v.IsStrictPositive()) {
    if (u.IsEven()) {
      shiftR(1, u.bits64);
      shiftL(1, s.bits64);
    } else if (v.IsEven()) {
      shiftR(1, v.bits64);
      shiftL(1, r.bits64);
    } else {
      x.Set(&u);
      x.Sub(&v);
      if (x.IsStrictPositive()) {
        shiftR(1, x.bits64);
        u.Set(&x);
        r.Add(&s);
        shiftL(1, s.bits64);
      } else {
        x.Neg();
        shiftR(1, x.bits64);
        v.Set(&x);
        s.Add(&r);
        shiftL(1, r.bits64);
      }
    }
    k++;
  }

  if (r.IsGreater(&Field_P))
    r.Sub(&Field_P);
  r.Neg();
  r.Add(&Field_P);

  for (int i = 0; i < k; i++) {
    if (r.IsEven()) {
      shiftR(1, r.bits64);
    } else {
      r.Add(&Field_P);
      shiftR(1, r.bits64);
    }
  }
  Set(&r);

#endif

#ifdef DRS62

  // Delayed right shift 62bits
  Int r;
  Int s;
  Int r0_P;
  Int s0_P;
  Int t1,t2,t3,t4;

  int64_t uu, uv, vu, vv;
  int64_t eta = -1;

  //printf("ModInv(%s)\n",GetBase16().c_str());


  //----------------------- First step (r,s) = (1,0)
  uu = 1; uv = 0;
  vu = 0; vv = 1;

  DivStep62(u.bits64[0],v.bits64[0],&eta,&uu,&uv,&vu,&vv);

  // Now update BigInt variables

  MatrixVecMul(u,v,uu,uv,vu,vv);
  LoadI64(t2,uv);
  LoadI64(t4,vv);

  // Compute multiple of P to add to s and r to make them multiple of 2^62
  uint64_t r0 = (t2.bits64[0] * MM64) & MSK62;
  uint64_t s0 = (t4.bits64[0] * MM64) & MSK62;
  r0_P.Mult(&Field_P,r0);
  s0_P.Mult(&Field_P,s0);
  r.Add(&t2,&r0_P);
  s.Add(&t4,&s0_P);

  // Right shift all variables by 62bits
  shiftR(62,u.bits64);
  shiftR(62,v.bits64);
  shiftR(62,r.bits64);
  shiftR(62,s.bits64);

  //----------------------- DivStep loop
  while (!v.IsZero()) {

    uu =  1; uv = 0;
    vu =  0; vv = 1;

    DivStep62(u.bits64[0],v.bits64[0],&eta,&uu,&uv,&vu,&vv);

    // Now update BigInt variables

    MatrixVecMul(u,v,uu,uv,vu,vv);
    MatrixVecMul(r,s,uu,uv,vu,vv);

    // Compute multiple of P to add to s and r to make them multiple of 2^62
    uint64_t r0 = (r.bits64[0] * MM64) & MSK62;
    uint64_t s0 = (s.bits64[0] * MM64) & MSK62;
    r0_P.Mult(&Field_P,r0);
    s0_P.Mult(&Field_P,s0);
    r.Add(&r0_P);
    s.Add(&s0_P);

    // Right shift all variables by 62bits
    shiftR(62, u.bits64);
    shiftR(62, v.bits64);
    shiftR(62, r.bits64);
    shiftR(62, s.bits64);

  }

  // u ends with -1 or 1
  if (u.IsNegative()) {
    // u = -1
    u.Neg();
    r.Neg();
  }
  if (!u.IsOne()) {
    // No inverse
    CLEAR();
    return;
  }

  while(r.IsNegative())
    r.Add(&Field_P);
  while(r.IsGreaterOrEqual(&Field_P))
    r.Sub(&Field_P);
  Set(&r);

#endif

}

// ------------------------------------------------

void Int::ModExp(Int *e) {

  Int base(this);
  SetInt32(1);
  uint32_t i = 0;

  uint32_t nbBit = e->GetBitLength();
  for(int i=0;i<(int)nbBit;i++) {
    if (e->GetBit(i))
      ModMul(&base);
    base.ModMul(&base);
  }

}

// ------------------------------------------------

void Int::ModMul(Int *a) {

  Int p;
  p.MontgomeryMult(a, this);
  MontgomeryMult(&_R2, &p);

}

// ------------------------------------------------

void Int::ModSquare(Int *a) {

  Int p;
  p.MontgomeryMult(a, a);
  MontgomeryMult(&_R2, &p);

}

// ------------------------------------------------

void Int::ModCube(Int *a) {

  Int p;
  Int p2;
  p.MontgomeryMult(a, a);
  p2.MontgomeryMult(&p, a);
  MontgomeryMult(&_R3, &p2);

}

// ------------------------------------------------

bool Int::HasSqrt() {

  // Euler's criterion
  Int e(&Field_P);
  Int a(this);
  e.SubOne();
  e.ShiftR(1);
  a.ModExp(&e);

  return a.IsOne();

}

// ------------------------------------------------

void Int::ModSqrt() {

  if (Field_P.IsEven()) {
    CLEAR();
    return;
  }

  if (!HasSqrt()) {
    CLEAR();
    return;
  }

  if ((Field_P.bits64[0] & 3) == 3) {

    Int e(&Field_P);
    e.AddOne();
    e.ShiftR(2);
    ModExp(&e);

  } else if ((Field_P.bits64[0] & 3) == 1) {

    int nbBit = Field_P.GetBitLength();

    // Tonelli Shanks
    uint64_t e=0;
    Int S(&Field_P);
    S.SubOne();
    while (S.IsEven()) {
      S.ShiftR(1);
      e++;
    }

    // Search smalest non-qresidue of P
    Int q((uint64_t)1);
    do {
      q.AddOne();
    }  while (q.HasSqrt());

    Int c(&q);
    c.ModExp(&S);

    Int t(this);
    t.ModExp(&S);

    Int r(this);
    Int ex(&S);
    ex.AddOne();
    ex.ShiftR(1);
    r.ModExp(&ex);

    uint64_t M = e;
    while (!t.IsOne()) {

      Int t2(&t);
      uint64_t i=0;
      while (!t2.IsOne()) {
        t2.ModSquare(&t2);
        i++;
      }

      Int b(&c);
      for(uint64_t j=0;j<M-i-1;j++)
        b.ModSquare(&b);
      M=i;
      c.ModSquare(&b);
      t.ModMul(&t,&c);
      r.ModMul(&r,&b);

    }

    Set(&r);

  }

}

// ------------------------------------------------

void Int::ModMul(Int *a, Int *b) {

  Int p;
  p.MontgomeryMult(a,b);
  MontgomeryMult(&_R2,&p);

}

// ------------------------------------------------

Int* Int::GetFieldCharacteristic() {
  return &Field_P;
}

// ------------------------------------------------

Int* Int::GetR() {
  return &_R;
}
Int* Int::GetR2() {
  return &_R2;
}
Int* Int::GetR3() {
  return &_R3;
}
Int* Int::GetR4() {
  return &_R4;
}

// ------------------------------------------------

void Int::SetupField(Int *n, Int *R, Int *R2, Int *R3, Int *R4) {

  // Size in number of 32bit word
  int nSize = n->GetSize();

  // Last digit inversions (Newton's iteration)
  {
    int64_t x, t;
    x = t = (int64_t)n->bits64[0];
    x = x * (2 - t * x);
    x = x * (2 - t * x);
    x = x * (2 - t * x);
    x = x * (2 - t * x);
    x = x * (2 - t * x);
    MM64 = (uint64_t)(-x);
    MM32 = (uint32_t)MM64;
  }
  Field_P.Set(n);

  // Size of Montgomery mult (64bits digit)
  Msize = nSize/2;

  // Compute few power of R
  // R = 2^(64*Msize) mod n
  Int Ri;
  Ri.MontgomeryMult(&_ONE, &_ONE); // Ri = R^-1
  _R.Set(&Ri);                     // R  = R^-1
  _R2.MontgomeryMult(&Ri, &_ONE);  // R2 = R^-2
  _R3.MontgomeryMult(&Ri, &Ri);    // R3 = R^-3
  _R4.MontgomeryMult(&_R3, &_ONE); // R4 = R^-4

  _R.ModInv();                     // R  = R
  _R2.ModInv();                    // R2 = R^2
  _R3.ModInv();                    // R3 = R^3
  _R4.ModInv();                    // R4 = R^4

  if (R)
    R->Set(&_R);

  if (R2)
    R2->Set(&_R2);

  if (R3)
    R3->Set(&_R3);

  if (R4)
    R4->Set(&_R4);

}
// ------------------------------------------------

uint64_t Int::AddC(Int *a) {

  unsigned char c = 0;
  c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
  c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
  c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
  c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
  c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);
#if NB64BLOCK > 5
  c = _addcarry_u64(c, bits64[5], a->bits64[5], bits64 + 5);
  c = _addcarry_u64(c, bits64[6], a->bits64[6], bits64 + 6);
  c = _addcarry_u64(c, bits64[7], a->bits64[7], bits64 + 7);
  c = _addcarry_u64(c, bits64[8], a->bits64[8], bits64 + 8);
#endif

  return c;

}

// ------------------------------------------------

void Int::AddAndShift(Int *a, Int *b, uint64_t cH) {

  unsigned char c = 0;
  c = _addcarry_u64(c, b->bits64[0], a->bits64[0], bits64 + 0);
  c = _addcarry_u64(c, b->bits64[1], a->bits64[1], bits64 + 0);
  c = _addcarry_u64(c, b->bits64[2], a->bits64[2], bits64 + 1);
  c = _addcarry_u64(c, b->bits64[3], a->bits64[3], bits64 + 2);
  c = _addcarry_u64(c, b->bits64[4], a->bits64[4], bits64 + 3);
#if NB64BLOCK > 5
  c = _addcarry_u64(c, b->bits64[5], a->bits64[5], bits64 + 4);
  c = _addcarry_u64(c, b->bits64[6], a->bits64[6], bits64 + 5);
  c = _addcarry_u64(c, b->bits64[7], a->bits64[7], bits64 + 6);
  c = _addcarry_u64(c, b->bits64[8], a->bits64[8], bits64 + 7);
#endif

  bits64[NB64BLOCK-1] = c + cH;

}

// ------------------------------------------------
void Int::MontgomeryMult(Int *a) {

  // Compute a*b*R^-1 (mod n),  R=2^k (mod n), k = Msize*64
  // a and b must be lower than n
  // See SetupField()

  Int t;
  Int pr;
  Int p;
  uint64_t ML;
  uint64_t c;

  // i = 0
  imm_umul(a->bits64, bits64[0], pr.bits64);
  ML = pr.bits64[0] * MM64;
  imm_umul(Field_P.bits64, ML, p.bits64);
  c = pr.AddC(&p);
  memcpy(t.bits64, pr.bits64 + 1, 8 * (NB64BLOCK - 1));
  t.bits64[NB64BLOCK - 1] = c;

  for (int i = 1; i < Msize; i++) {

    imm_umul(a->bits64, bits64[i], pr.bits64);
    ML = (pr.bits64[0] + t.bits64[0]) * MM64;
    imm_umul(Field_P.bits64, ML, p.bits64);
	  c = pr.AddC(&p);
    t.AddAndShift(&t, &pr, c);

  }

  p.Sub(&t,&Field_P);
  if (p.IsPositive())
    Set(&p);
  else
    Set(&t);

}

void Int::MontgomeryMult(Int *a, Int *b) {

  // Compute a*b*R^-1 (mod n),  R=2^k (mod n), k = Msize*64
  // a and b must be lower than n
  // See SetupField()

  Int pr;
  Int p;
  uint64_t ML;
  uint64_t c;

  // i = 0
  imm_umul(a->bits64, b->bits64[0], pr.bits64);
  ML = pr.bits64[0] * MM64;
  imm_umul(Field_P.bits64, ML, p.bits64);
  c = pr.AddC(&p);
  memcpy(bits64,pr.bits64 + 1,8*(NB64BLOCK-1));
  bits64[NB64BLOCK-1] = c;

  for (int i = 1; i < Msize; i++) {

    imm_umul(a->bits64, b->bits64[i], pr.bits64);
    ML = (pr.bits64[0] + bits64[0]) * MM64;
    imm_umul(Field_P.bits64, ML, p.bits64);
	  c = pr.AddC(&p);
    AddAndShift(this, &pr, c);

  }

  p.Sub(this, &Field_P);
  if (p.IsPositive())
    Set(&p);

}


// SecpK1 specific section -----------------------------------------------------------------------------

void Int::ModMulK1(Int *a, Int *b) {
/*
#ifndef WIN64
#if (__GNUC__ > 7) || (__GNUC__ == 7 && (__GNUC_MINOR__ > 2))
  unsigned char c;
#else
  #warning "GCC lass than 7.3 detected, upgrade gcc to get best perfromance"
  volatile unsigned char c;
#endif
#else
  unsigned char c;
#endif
*/
  unsigned char c;

  uint64_t ah, al;
  uint64_t t[5];
  uint64_t r512[8];
  r512[5] = 0;
  r512[6] = 0;
  r512[7] = 0;

  // 256*256 multiplier
  imm_umul(a->bits64, b->bits64[0], r512);
  imm_umul(a->bits64, b->bits64[1], t);
  c = _addcarry_u64(0, r512[1], t[0], r512 + 1);
  c = _addcarry_u64(c, r512[2], t[1], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[2], r512 + 3);
  c = _addcarry_u64(c, r512[4], t[3], r512 + 4);
  c = _addcarry_u64(c, r512[5], t[4], r512 + 5);
  imm_umul(a->bits64, b->bits64[2], t);
  c = _addcarry_u64(0, r512[2], t[0], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[1], r512 + 3);
  c = _addcarry_u64(c, r512[4], t[2], r512 + 4);
  c = _addcarry_u64(c, r512[5], t[3], r512 + 5);
  c = _addcarry_u64(c, r512[6], t[4], r512 + 6);
  imm_umul(a->bits64, b->bits64[3], t);
  c = _addcarry_u64(0, r512[3], t[0], r512 + 3);
  c = _addcarry_u64(c, r512[4], t[1], r512 + 4);
  c = _addcarry_u64(c, r512[5], t[2], r512 + 5);
  c = _addcarry_u64(c, r512[6], t[3], r512 + 6);
  c = _addcarry_u64(c, r512[7], t[4], r512 + 7);

  // Reduce from 512 to 320
  imm_umul(r512 + 4, 0x1000003D1ULL, t);
  c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
  c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
  c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

  // Reduce from 320 to 256
  // No overflow possible here t[4]+c<=0x1000003D1ULL
  al = _umul128(t[4] + c, 0x1000003D1ULL, &ah);
  c = _addcarry_u64(0, r512[0], al, bits64 + 0);
  c = _addcarry_u64(c, r512[1], ah, bits64 + 1);
  c = _addcarry_u64(c, r512[2], 0ULL, bits64 + 2);
  c = _addcarry_u64(c, r512[3], 0ULL, bits64 + 3);

  // Probability of carry here or that this>P is very very unlikely
  bits64[4] = 0;

}

void Int::ModMulK1(Int *a) {
/*
#ifndef WIN64
#if (__GNUC__ > 7) || (__GNUC__ == 7 && (__GNUC_MINOR__ > 2))
  unsigned char c;
#else
  #warning "GCC lass than 7.3 detected, upgrade gcc to get best perfromance"
  volatile unsigned char c;
#endif
#else
  unsigned char c;
#endif
*/
  unsigned char c;

  uint64_t ah, al;
  uint64_t t[5];
  uint64_t r512[8];
  r512[5] = 0;
  r512[6] = 0;
  r512[7] = 0;

  // 256*256 multiplier
  imm_umul(a->bits64, bits64[0], r512);
  imm_umul(a->bits64, bits64[1], t);
  c = _addcarry_u64(0, r512[1], t[0], r512 + 1);
  c = _addcarry_u64(c, r512[2], t[1], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[2], r512 + 3);
  c = _addcarry_u64(c, r512[4], t[3], r512 + 4);
  c = _addcarry_u64(c, r512[5], t[4], r512 + 5);
  imm_umul(a->bits64, bits64[2], t);
  c = _addcarry_u64(0, r512[2], t[0], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[1], r512 + 3);
  c = _addcarry_u64(c, r512[4], t[2], r512 + 4);
  c = _addcarry_u64(c, r512[5], t[3], r512 + 5);
  c = _addcarry_u64(c, r512[6], t[4], r512 + 6);
  imm_umul(a->bits64, bits64[3], t);
  c = _addcarry_u64(0, r512[3], t[0], r512 + 3);
  c = _addcarry_u64(c, r512[4], t[1], r512 + 4);
  c = _addcarry_u64(c, r512[5], t[2], r512 + 5);
  c = _addcarry_u64(c, r512[6], t[3], r512 + 6);
  c = _addcarry_u64(c, r512[7], t[4], r512 + 7);

  // Reduce from 512 to 320
  imm_umul(r512 + 4, 0x1000003D1ULL, t);
  c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
  c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
  c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

  // Reduce from 320 to 256
  // No overflow possible here t[4]+c<=0x1000003D1ULL
  al = _umul128(t[4] + c, 0x1000003D1ULL, &ah);
  c = _addcarry_u64(0, r512[0], al, bits64 + 0);
  c = _addcarry_u64(c, r512[1], ah, bits64 + 1);
  c = _addcarry_u64(c, r512[2], 0, bits64 + 2);
  c = _addcarry_u64(c, r512[3], 0, bits64 + 3);
  // Probability of carry here or that this>P is very very unlikely
  bits64[4] = 0;

}

void Int::ModSquareK1(Int *a) {
/*
#ifndef WIN64
#if (__GNUC__ > 7) || (__GNUC__ == 7 && (__GNUC_MINOR__ > 2))
  unsigned char c;
#else
  #warning "GCC lass than 7.3 detected, upgrade gcc to get best perfromance"
  volatile unsigned char c;
#endif
#else
  unsigned char c;
#endif
*/
  unsigned char c;

  uint64_t r512[8];
  uint64_t u10, u11;
  uint64_t t1;
  uint64_t t2;
  uint64_t t[5];


  //k=0
  r512[0] = _umul128(a->bits64[0], a->bits64[0], &t[1]);

  //k=1
  t[3] = _umul128(a->bits64[0], a->bits64[1], &t[4]);
  c = _addcarry_u64(0, t[3], t[3], &t[3]);
  c = _addcarry_u64(c, t[4], t[4], &t[4]);
  c = _addcarry_u64(c,  0,  0, &t1);
  c = _addcarry_u64(0, t[1], t[3], &t[3]);
  c = _addcarry_u64(c, t[4],  0, &t[4]);
  c = _addcarry_u64(c, t1,  0, &t1);
  r512[1] = t[3];

  //k=2
  t[0] = _umul128(a->bits64[0], a->bits64[2], &t[1]);
  c = _addcarry_u64(0, t[0], t[0], &t[0]);
  c = _addcarry_u64(c, t[1], t[1], &t[1]);
  c = _addcarry_u64(c,  0,  0, &t2);

  u10 = _umul128(a->bits64[1], a->bits64[1], &u11);
  c = _addcarry_u64(0, t[0] , u10, &t[0]);
  c = _addcarry_u64(c, t[1] , u11, &t[1]);
  c = _addcarry_u64(c, t2 ,   0, &t2);
  c = _addcarry_u64(0, t[0], t[4], &t[0]);
  c = _addcarry_u64(c, t[1], t1, &t[1]);
  c = _addcarry_u64(c, t2, 0, &t2);
  r512[2] = t[0];

  //k=3
  t[3] = _umul128(a->bits64[0], a->bits64[3], &t[4]);
  u10 = _umul128(a->bits64[1], a->bits64[2], &u11);

  c = _addcarry_u64(0, t[3], u10, &t[3]);
  c = _addcarry_u64(c, t[4], u11, &t[4]);
  c = _addcarry_u64(c,  0,   0, &t1);
  t1 += t1;
  c = _addcarry_u64(0, t[3], t[3], &t[3]);
  c = _addcarry_u64(c, t[4], t[4], &t[4]);
  c = _addcarry_u64(c, t1, 0, &t1);
  c = _addcarry_u64(0, t[3], t[1], &t[3]);
  c = _addcarry_u64(c, t[4], t2, &t[4]);
  c = _addcarry_u64(c, t1, 0, &t1);
  r512[3] = t[3];

  //k=4
  t[0] = _umul128(a->bits64[1], a->bits64[3], &t[1]);
  c = _addcarry_u64(0, t[0], t[0], &t[0]);
  c = _addcarry_u64(c, t[1], t[1], &t[1]);
  c = _addcarry_u64(c, 0, 0, &t2);

  u10 = _umul128(a->bits64[2], a->bits64[2], &u11);
  c = _addcarry_u64(0, t[0], u10, &t[0]);
  c = _addcarry_u64(c, t[1], u11, &t[1]);
  c = _addcarry_u64(c, t2, 0, &t2);
  c = _addcarry_u64(0, t[0], t[4], &t[0]);
  c = _addcarry_u64(c, t[1], t1, &t[1]);
  c = _addcarry_u64(c, t2,  0, &t2);
  r512[4] = t[0];

  //k=5
  t[3] = _umul128(a->bits64[2], a->bits64[3], &t[4]);
  c = _addcarry_u64(0, t[3], t[3], &t[3]);
  c = _addcarry_u64(c, t[4], t[4], &t[4]);
  c = _addcarry_u64(c, 0, 0, &t1);
  c = _addcarry_u64(0, t[3], t[1], &t[3]);
  c = _addcarry_u64(c, t[4], t2, &t[4]);
  c = _addcarry_u64(c, t1,  0, &t1);
  r512[5] = t[3];

  //k=6
  t[0] = _umul128(a->bits64[3], a->bits64[3], &t[1]);
  c = _addcarry_u64(0, t[0], t[4], &t[0]);
  c = _addcarry_u64(c, t[1], t1, &t[1]);
  r512[6] = t[0];

  //k=7
  r512[7] = t[1];

  // Reduce from 512 to 320
  imm_umul(r512 + 4, 0x1000003D1ULL, t);
  c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
  c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
  c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
  c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

  // Reduce from 320 to 256
  // No overflow possible here t[4]+c<=0x1000003D1ULL
  u10 = _umul128(t[4] + c, 0x1000003D1ULL, &u11);
  c = _addcarry_u64(0, r512[0], u10, bits64 + 0);
  c = _addcarry_u64(c, r512[1], u11, bits64 + 1);
  c = _addcarry_u64(c, r512[2], 0, bits64 + 2);
  c = _addcarry_u64(c, r512[3], 0, bits64 + 3);
  // Probability of carry here or that this>P is very very unlikely
  bits64[4] = 0;

}

static Int _R2o;                               // R^2 for SecpK1 order modular mult
static uint64_t MM64o = 0x4B0DFF665588B13FULL; // 64bits lsb negative inverse of SecpK1 order
static Int *_O;                                // SecpK1 order

void Int::InitK1(Int *order) {
  _O = order;
  _R2o.SetBase16("9D671CD581C69BC5E697F5E45BCD07C6741496C20E7CF878896CF21467D7D140");
}

void Int::ModAddK1order(Int *a, Int *b) {
  Add(a,b);
  Sub(_O);
  if (IsNegative())
    Add(_O);
}

void Int::ModAddK1order(Int *a) {
  Add(a);
  Sub(_O);
  if(IsNegative())
    Add(_O);
}

void Int::ModSubK1order(Int *a) {
  Sub(a);
  if(IsNegative())
    Add(_O);
}

void Int::ModNegK1order() {
  Neg();
  Add(_O);
}

uint32_t Int::ModPositiveK1() {

  Int N(this);
  Int D(this);
  N.ModNeg();
  D.Sub(&N);
  if(D.IsNegative()) {
    return 0;
  } else {
    Set(&N);
    return 1;
  }

}


void Int::ModMulK1order(Int *a) {

  Int t;
  Int pr;
  Int p;
  uint64_t ML;
  uint64_t c;

  imm_umul(a->bits64, bits64[0], pr.bits64);
  ML = pr.bits64[0] * MM64o;
  imm_umul(_O->bits64, ML, p.bits64);
  c = pr.AddC(&p);
  memcpy(t.bits64, pr.bits64 + 1, 32);
  t.bits64[4] = c;

  for (int i = 1; i < 4; i++) {

    imm_umul(a->bits64, bits64[i], pr.bits64);
    ML = (pr.bits64[0] + t.bits64[0]) * MM64o;
    imm_umul(_O->bits64, ML, p.bits64);
    c = pr.AddC(&p);
    t.AddAndShift(&t, &pr, c);

  }

  p.Sub(&t, _O);
  if (p.IsPositive())
    Set(&p);
  else
    Set(&t);


  // Normalize

  imm_umul(_R2o.bits64, bits64[0], pr.bits64);
  ML = pr.bits64[0] * MM64o;
  imm_umul(_O->bits64, ML, p.bits64);
  c = pr.AddC(&p);
  memcpy(t.bits64, pr.bits64 + 1, 32);
  t.bits64[4] = c;

  for (int i = 1; i < 4; i++) {

    imm_umul(_R2o.bits64, bits64[i], pr.bits64);
    ML = (pr.bits64[0] + t.bits64[0]) * MM64o;
    imm_umul(_O->bits64, ML, p.bits64);
    c = pr.AddC(&p);
    t.AddAndShift(&t, &pr, c);

  }

  p.Sub(&t, _O);
  if (p.IsPositive())
    Set(&p);
  else
    Set(&t);

}
