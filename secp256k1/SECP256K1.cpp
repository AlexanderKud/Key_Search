#include "SECP256k1.h"
#include <string.h>
#include <cstdint>
#include <iostream>
#include <gmp.h>

Secp256K1::Secp256K1() {
}

void Secp256K1::Init() {

  // Prime for the finite field
  //Int P;
  P.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");
  // Set up field
  Int::SetupField(&P);
  // Generator point and order
  G.x.SetBase16("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798");
  G.y.SetBase16("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8");
  G.z.SetBase16("0000000000000000000000000000000000000000000000000000000000000001");
  order.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141");
  Int::InitK1(&order);
  // Compute Generator table
  Point N(G);
  for(int i = 0; i < 32; i++) {
    GTable[i * 256] = N;
    N = DoublePoint(N);
    for (int j = 1; j < 255; j++) {
      GTable[i * 256 + j] = N;
      N = AddPoints(N, GTable[i * 256]);
    }
  }

}

Point Secp256K1::ScalarMultiplication(Int *privKey) {

  int i = 0;
  uint8_t b;
  Point Q;
  Q.Clear();

  // Search first significant byte
  for (i = 0; i < 32; i++) {
    b = privKey->GetByte(i);
    if(b) break;
  }
  Q = GTable[256 * i + (b-1)];
  i++;

  for(; i < 32; i++) {
    b = privKey->GetByte(i);
    if(b) Q = Add2(Q, GTable[256 * i + (b-1)]);
  }

  Q.Reduce();
  return Q;

}

Point Secp256K1::PointMultiplication(Point &P, Int *scalar) {
  Point R, T;
  int  no_of_bits, loop;
  no_of_bits = scalar->GetBitLength();
  R.Set(P); R.z.SetInt32(1);
  T.Set(P); T.z.SetInt32(1);
  for(loop = no_of_bits - 2; loop >= 0; loop--) {
      R = Double(R);
      if(scalar->GetBit(loop)) { R = Add2(R, T); }        
  }
  R.Reduce();
  return R;
}

Point Secp256K1::PointDivision(Point &P, Int *scalar) {
  Point A;
  Int mod_inv;
  mpz_t N, d, s;
  mpz_inits(N, d, s, NULL);
  mpz_set_str(N, order.GetBase10().c_str(), 0);
  mpz_set_str(s, scalar->GetBase10().c_str(), 0);
  mpz_invert(d, s, N);
  mod_inv.SetBase10(mpz_get_str(NULL, 10, d));
  mpz_clears(N, d, s, NULL);
  A = PointMultiplication(P, &mod_inv);
  return A;
}

uint8_t Secp256K1::GetByte(std::string &str, int idx) {

  char tmp[3];
  int  val;

  tmp[0] = str.data()[2 * idx];
  tmp[1] = str.data()[2 * idx + 1];
  tmp[2] = 0;

  sscanf(tmp, "%X", &val);
  return (uint8_t)val;

}

Point Secp256K1::ParsePublicKeyHex(std::string str) {

  Point ret;
  ret.Clear();

  uint8_t type = GetByte(str, 0);

  switch (type) {

    case 0x02:
      for (int i = 0; i < 32; i++)
        ret.x.SetByte(31 - i, GetByte(str, i + 1));
      ret.y = GetY(ret.x, true);
      break;

    case 0x03:
      for (int i = 0; i < 32; i++)
        ret.x.SetByte(31 - i, GetByte(str, i + 1));
      ret.y = GetY(ret.x, false);
      break;

  }

  ret.z.SetInt32(1);
  return ret;

}

std::string Secp256K1::GetPublicKeyHex(Point &pubKey) {

  unsigned char publicKeyBytes[128];
  char tmp[3];
  std::string ret;
  publicKeyBytes[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
  pubKey.x.Get32Bytes(publicKeyBytes + 1);
  for (int i = 0; i < 33; i++) {
    sprintf(tmp, "%02x", (int)publicKeyBytes[i]);
    ret.append(tmp);
  }

  return ret;
}

std::string Secp256K1::GetXHex(Int* x, int length) {
  unsigned char publicKeyBytes[33];
  char tmp[3];
  std::string ret;
  x->Get32Bytes(publicKeyBytes);
  for (int i = 1; i < length; i++) {
    sprintf(tmp, "%02x", (int)publicKeyBytes[i]);
    ret.append(tmp);
  }

  return ret;
}

Point Secp256K1::AddPoints(Point &p1,Point &p2) {

  Int _s, dx, dy;
  Point r;
  r.z.SetInt32(1);

  dy.ModSub(&p2.y, &p1.y);
  dx.ModSub(&p2.x, &p1.x);
  dx.ModInv();
  _s.ModMulK1(&dy, &dx);     // s = (p2.y-p1.y)*inverse(p2.x-p1.x);

  r.x.ModSquareK1(&_s);       // _p = pow2(s)
  r.x.ModSub(&p1.x);
  r.x.ModSub(&p2.x);       // rx = pow2(s) - p1.x - p2.x;

  r.y.ModSub(&p2.x, &r.x);
  r.y.ModMulK1(&_s);
  r.y.ModSub(&p2.y);       // ry = - p2.y - s*(ret.x-p2.x);

  return r;

}

Point Secp256K1::Add2(Point &p1, Point &p2) {
    Int u;
    Int v;
    Int u1;
    Int v1;
    Int vs2;
    Int vs3;
    Int us2;
    Int a;
    Int us2w;
    Int vs2v2;
    Int vs3u2;
    Int _2vs2v2;
    Point r;

    u1.ModMulK1(&p2.y, &p1.z);
    v1.ModMulK1(&p2.x, &p1.z);
    u.ModSub(&u1, &p1.y);
    v.ModSub(&v1, &p1.x);
    us2.ModSquareK1(&u);
    vs2.ModSquareK1(&v);
    vs3.ModMulK1(&vs2, &v);
    us2w.ModMulK1(&us2, &p1.z);
    vs2v2.ModMulK1(&vs2, &p1.x);
    _2vs2v2.ModAdd(&vs2v2, &vs2v2);
    a.ModSub(&us2w, &vs3);
    a.ModSub(&_2vs2v2);

    r.x.ModMulK1(&v, &a);

    vs3u2.ModMulK1(&vs3, &p1.y);
    r.y.ModSub(&vs2v2, &a);
    r.y.ModMulK1(&r.y, &u);
    r.y.ModSub(&vs3u2);

    r.z.ModMulK1(&vs3, &p1.z);

    return r;
}

Point Secp256K1::Add(Point &p1, Point &p2) {

  Int u;
  Int v;
  Int u1;
  Int u2;
  Int v1;
  Int v2;
  Int vs2;
  Int vs3;
  Int us2;
  Int w;
  Int a;
  Int us2w;
  Int vs2v2;
  Int vs3u2;
  Int _2vs2v2;
  Int x3;
  Int vs3y1;
  Point r;

  /*
  U1 = Y2 * Z1
  U2 = Y1 * Z2
  V1 = X2 * Z1
  V2 = X1 * Z2
  if (V1 == V2)
    if (U1 != U2)
      return POINT_AT_INFINITY
    else
      return POINT_DOUBLE(X1, Y1, Z1)
  U = U1 - U2
  V = V1 - V2
  W = Z1 * Z2
  A = U ^ 2 * W - V ^ 3 - 2 * V ^ 2 * V2
  X3 = V * A
  Y3 = U * (V ^ 2 * V2 - A) - V ^ 3 * U2
  Z3 = V ^ 3 * W
  return (X3, Y3, Z3)
  */

  u1.ModMulK1(&p2.y,&p1.z);
  u2.ModMulK1(&p1.y,&p2.z);
  v1.ModMulK1(&p2.x,&p1.z);
  v2.ModMulK1(&p1.x,&p2.z);
  u.ModSub(&u1,&u2);
  v.ModSub(&v1,&v2);
  w.ModMulK1(&p1.z,&p2.z);
  us2.ModSquareK1(&u);
  vs2.ModSquareK1(&v);
  vs3.ModMulK1(&vs2,&v);
  us2w.ModMulK1(&us2,&w);
  vs2v2.ModMulK1(&vs2,&v2);
  _2vs2v2.ModAdd(&vs2v2,&vs2v2);
  a.ModSub(&us2w,&vs3);
  a.ModSub(&_2vs2v2);

  r.x.ModMulK1(&v,&a);

  vs3u2.ModMulK1(&vs3,&u2);
  r.y.ModSub(&vs2v2,&a);
  r.y.ModMulK1(&r.y,&u);
  r.y.ModSub(&vs3u2);

  r.z.ModMulK1(&vs3,&w);

  return r;
}

Point Secp256K1::DoublePoint(Point &p) {

  Int _s;
  Int _p;
  Int a;
  Point r;
  r.z.SetInt32(1);

  _s.ModMulK1(&p.x,&p.x);
  _p.ModAdd(&_s,&_s);
  _p.ModAdd(&_s);

  a.ModAdd(&p.y,&p.y);
  a.ModInv();
  _s.ModMulK1(&_p,&a);     // s = (3*pow2(p.x))*inverse(2*p.y);

  _p.ModMulK1(&_s,&_s);
  a.ModAdd(&p.x,&p.x);
  a.ModNeg();
  r.x.ModAdd(&a,&_p);    // rx = pow2(s) + neg(2*p.x);

  a.ModSub(&r.x,&p.x);

  _p.ModMulK1(&a,&_s);
  r.y.ModAdd(&_p,&p.y);
  r.y.ModNeg();           // ry = neg(p.y + s*(ret.x+neg(p.x)));

  return r;
}

Point Secp256K1::Double(Point &p) {


  /*
  if (Y == 0)
    return POINT_AT_INFINITY
    W = a * Z ^ 2 + 3 * X ^ 2
    S = Y * Z
    B = X * Y*S
    H = W ^ 2 - 8 * B
    X' = 2*H*S
    Y' = W*(4*B - H) - 8*Y^2*S^2
    Z' = 8*S^3
    return (X', Y', Z')
  */

  Int z2;
  Int x2;
  Int _3x2;
  Int w;
  Int s;
  Int s2;
  Int b;
  Int _8b;
  Int _8y2s2;
  Int y2;
  Int h;
  Point r;

  z2.ModSquareK1(&p.z);
  z2.SetInt32(0); // a=0
  x2.ModSquareK1(&p.x);
  _3x2.ModAdd(&x2,&x2);
  _3x2.ModAdd(&x2);
  w.ModAdd(&z2,&_3x2);
  s.ModMulK1(&p.y,&p.z);
  b.ModMulK1(&p.y,&s);
  b.ModMulK1(&p.x);
  h.ModSquareK1(&w);
  _8b.ModAdd(&b,&b);
  _8b.ModDouble();
  _8b.ModDouble();
  h.ModSub(&_8b);

  r.x.ModMulK1(&h,&s);
  r.x.ModAdd(&r.x);

  s2.ModSquareK1(&s);
  y2.ModSquareK1(&p.y);
  _8y2s2.ModMulK1(&y2,&s2);
  _8y2s2.ModDouble();
  _8y2s2.ModDouble();
  _8y2s2.ModDouble();

  r.y.ModAdd(&b,&b);
  r.y.ModAdd(&r.y,&r.y);
  r.y.ModSub(&h);
  r.y.ModMulK1(&w);
  r.y.ModSub(&_8y2s2);

  r.z.ModMulK1(&s2,&s);
  r.z.ModDouble();
  r.z.ModDouble();
  r.z.ModDouble();

  return r;
}

Point Secp256K1::SubtractPoints(Point &p1, Point &p2) {
  Point Q1, Q2;
  Q1.Set(p2);
  Q1.y.ModNeg();
  Q1.z.SetInt32(1);
  Q2 = AddPoints(p1, Q1);
  return Q2;
}

Int Secp256K1::GetY(Int x, bool isEven) {

  Int _s;
  Int _p;

  _s.ModSquareK1(&x);
  _p.ModMulK1(&_s,&x);
  _p.ModAdd(7);
  _p.ModSqrt();

  if(!_p.IsEven() && isEven) {
    _p.ModNeg();
  }
  else if(_p.IsEven() && !isEven) {
    _p.ModNeg();
  }

  return _p;

}

bool Secp256K1::EC(Point &p) {

  Int _s;
  Int _p;

  _s.ModSquareK1(&p.x);
  _p.ModMulK1(&_s,&p.x);
  _p.ModAdd(7);
  _s.ModMulK1(&p.y,&p.y);
  _s.ModSub(&_p);

  return _s.IsZero(); // ( ((pow2(y) - (pow3(x) + 7)) % P) == 0 );

}
