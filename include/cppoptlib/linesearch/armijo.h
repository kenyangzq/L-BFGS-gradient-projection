// CppNumericalSolver
#ifndef ARMIJO_H_
#define ARMIJO_H_

#include "../meta.h"

namespace cppoptlib {

template<typename T, typename P, int Ord>
class Armijo {

 public:

    /**
     * @brief use Armijo Rule for (weak) Wolfe conditiions
     * @details [long description]
     *
     * @param searchDir search direction for next update step
     * @param objFunc handle to problem
     *
     * @return step-width
     */
public:
    
    static void linesearch(Vector<T> & x0, const Vector<T> & searchDir, const P &objFunc, const T alpha_init = 1) {
        
        Vector<T> x = x0;
        
        // evaluate phi(0)
        T phi0 = objFunc.value(x0);
        
        
        
        // evaluate phi'(0)
        Vector<T> grad(x.rows());
        objFunc.gradient(x, grad);
        
        
        T phi0_dash = searchDir.dot(grad);
        
        T alpha = alpha_init;
        T phi = 0;
        Vector<T> x_candidate (x.rows());
        Vector<T> x_min (x.rows());
        Vector<T> grad2(x.rows());
        
        double min = phi0;
        
        
        // 200 guesses
        for(size_t iter = 0; iter < 10; ++iter) {
            
            // new guess for phi(alpha)
            x_candidate = x + alpha * searchDir;
            objFunc.pullBack(x_candidate);
            phi = objFunc.value(x_candidate);
            
            
            // decrease condition invalid --> shrink interval
            if (phi > phi0 + 0.0001 * alpha * phi0_dash) {
                alpha *= 0.3;
                
            } else {
                
                // valid decrease --> test strong wolfe condition
                
                objFunc.gradient(x_candidate, grad2);
                const T phi_dash = searchDir.dot(grad2);
                
                // curvature condition invalid ?
                if (phi_dash < 0.9 * phi0_dash) {
                    // increase interval
                    alpha *= 4.5;
                } else if (phi < phi0){
                    // both condition are valid --> we are happy
                    x0 = x_candidate;
                    printf("phi-dash: %f\n", phi_dash);
                    return;
                }
                
                if (phi < min) {
                    min = phi;
                    x_min = x_candidate;
                }
                
                
                
            }
        }
        if (min < phi0) {
            x0 = x_min;
        }
        return;
    }


};

//template<typename T, typename P>
//class Armijo<T, P, 2> {
//
// public:
//
//    /**
//     * @brief use Armijo Rule for (weak) Wolfe conditiions
//     * @details [long description]
//     *
//     * @param searchDir search direction for next update step
//     * @param objFunc handle to problem
//     *
//     * @return step-width
//     */
//    static T linesearch(const Vector<T> & x, const Vector<T> & searchDir, P &objFunc) {
//        const T c = 0.2;
//        const T rho = 0.9;
//        T alpha = 1.0;
//
//        T f = objFunc.value(x + alpha * searchDir);
//        const T f_in = objFunc.value(x);
//        const Matrix<T>  hessian(x.rows(), x.rows());
//        objFunc.hessian(x, hessian);
//        Vector<T> grad(x.rows());
//        objFunc.gradient(x, grad);
//        const T Cache = c * grad.dot(searchDir) + 0.5 * c*c * searchDir.transpose() * (hessian * searchDir);
//
//        while(f > f_in + alpha * Cache) {
//            alpha *= rho;
//            f = objFunc.value(x + alpha * searchDir);
//        }
//        return alpha;
//    }
//
//};

}

#endif /* ARMIJO_H_ */
