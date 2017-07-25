#ifndef vec3_hpp
#define vec3_hpp

#include <iostream>
#include <vector>

class vector3{
    public:        
        vector3(float x, float y, float z);
        vector3();
        float x(void);
        float y(void);
        float z(void);
        void normalize(void);
        void setValue(float x, float y, float z);
        static inline float dotproduct(vector3 vec1, vector3 vec2){
            return vec1.x()*vec2.x()+vec1.y()*vec2.y()+vec1.z()*vec2.z();
        }
       static inline vector3 crossproduct(vector3 u, vector3 v){ 
        float i,j,k;
            i = u.y()*v.z() - u.z()*v.y();
            j = u.z()*v.x() - u.x()*v.z();
            k = u.x()*v.y() - u.y()*v.x();
            vector3 vec(i, j, k);
            return vec;
        }
        static vector3 vec_scal_mult(float c, vector3 v){
            float x,y,z;
            x = c*v.x();
            y=c*v.y();
            z=c*v.z();
            vector3 vec(x,y,z);
            return vec;
        }
        static inline vector3 add(vector3 v1, vector3 v2){
            float x,y,z;
            x = v1.x()+v2.x();
            y = v1.y()+v2.y();
            z = v1.z()+v2.z();
            vector3 vec(x,y,z);
            return vec;
        }
          static inline vector3 subtract(vector3 v1, vector3 v2){
            float x,y,z;
            x = v1.x()-v2.x();
            y = v1.y()-v2.y();
            z = v1.z()-v2.z();
            vector3 vec(x,y,z);
            return vec;
        }
        static inline vector3 add3(vector3 v1, vector3 v2, vector3 v3){
            float x,y,z;
            x = v1.x()+v2.x()+v3.x();
            y = v1.y()+v2.y()+v3.y();
            z = v1.z()+v2.z()+v3.z();
            vector3 vec(x,y,z);
            return vec;
        }
    private:    
        float x_val;
        float y_val;
        float z_val;
};
#endif