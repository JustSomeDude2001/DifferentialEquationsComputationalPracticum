class DifferentialEquation
{
private:    
    float x_0;
    float y_0;
    float (*y_exact)(float);
    float (*y_sharp)(float, float);
public:
    
    float getX_0() {
        return x_0;
    }

    float getY_0() {
        return y_0;
    }

    (float*)(float, float) getY_exact() {
    }

    DifferentialEquation();
    ~DifferentialEquation();
};
