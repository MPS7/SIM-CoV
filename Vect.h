
#ifndef VECT_H
#define VECT_H

#ifndef WWWW
#define WWWW
#define WWWWW(a) {if(((INOUT)?fread(&(a),sizeof(a),1,file):fwrite(&(a),sizeof(a),1,file))==0)return false;}
#endif //WWWW
;
class Vect2D
{
 public:
  Vect2D(){x=0; y=0;}
  Vect2D(double x_, double y_){x=x_; y=y_;}
  ~Vect2D(){};

  void operator= (Vect2D v) {x=v.x; y=v.y;}
  void operator= (double z) {x=z; y=z;}
  Vect2D operator+ (Vect2D v) {return Vect2D(x+v.x,y+v.y);}
  void operator+= (Vect2D v) {x+=v.x; y+=v.y;}
  Vect2D operator- (Vect2D v) {return Vect2D(x-v.x,y-v.y);}
  void operator-= (Vect2D v) {x-=v.x; y-=v.y;}
  Vect2D operator* (double lambda) {return Vect2D(lambda*x,lambda*y);}
  Vect2D operator/ (double lambda) {return Vect2D(x/lambda,y/lambda);}


  double x;
  double y;
};

#endif //VECT_H
