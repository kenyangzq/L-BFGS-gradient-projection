// CppNumericalSolver
#ifndef WOLFERULE_H_
#define WOLFERULE_H_

#include "../meta.h"

namespace cppoptlib {
    
    /**
     * @brief this tries to guess the correct stepwith in a bisection-style
     * @details WARNING: THIS IS QUITE HACKY. TEST ARMIJO before.
     *
     * @tparam T scalar type
     * @tparam P problem type
     * @tparam Ord order of solver
     */
    template<typename T, typename P, int Ord>
    class WolfeHeuristic {
        
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
            
            // 200 guesses
            for(size_t iter = 0; iter < 20; ++iter) {
                printf("alpha: %f\t", alpha);
                
                // new guess for phi(alpha)
                x_candidate = x + alpha * searchDir;
                objFunc.pullBack(x_candidate);
                phi = objFunc.value(x_candidate);
            
                // decrease condition invalid --> shrink interval
                if (phi > phi0 + 0.0001 * alpha * phi0_dash) {
                    alpha *= 0.4;
                    
                } else {
    
                    // valid decrease --> test strong wolfe condition
                    objFunc.gradient(x_candidate, grad2);
                    const T phi_dash = searchDir.dot(grad2);
    
                    // curvature condition invalid ?
                    if (phi_dash < 0.9 * phi0_dash) {
                        // increase interval
                        alpha *= 3; 
                    } else if (phi < phi0){
                        // both condition are valid --> we are happy
                        x0 = x_candidate;
                        printf("phi-dash: %f\n", phi_dash);
                        return;
                    }
                }
            }
            if (phi < phi0) {
                x0 = x_candidate;
            }
            return;
        }
    };
}

    
    
#endif /* WOLFERULE_H_ */
