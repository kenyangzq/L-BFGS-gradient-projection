// CppNumericalSolver
#include <iostream>
#include <Eigen/LU>
#include <iomanip>
#include <fstream>
#include "isolver.h"
#include "../linesearch/armijo.h"
#include "../linesearch/morethuente.h"
#include "../linesearch/wolfeheuristic.h"

#ifndef LBFGSSOLVER_H_
#define LBFGSSOLVER_H_

namespace std {
    void writeFiles (ofstream & outputfile, string name, cppoptlib::Vector<double> V, int dim){
        
        int c = dim-1;
        
        outputfile << setprecision(6);
        outputfile << fixed;
        
        cppoptlib::Vector<double> tmp2(dim);
        cppoptlib::Vector<double> tmp(c);
        
        for (int i =0; i < V.rows()/c; i++) {
            tmp = V.segment(i*c, c);

            tmp2(0) = cos(tmp(0)) * sin(tmp(1));
            tmp2(1) = sin(tmp(0)) * sin(tmp(1));
            tmp2(2) = cos(tmp(1));
            
            outputfile << tmp2(0) << "\t" << tmp2(1) << "\t" << tmp2(2) << "\n";
        }
    }
}







namespace cppoptlib {

template<typename T>
class LbfgsSolver : public ISolver<T, 1> {
  public:
    
    size_t numFile = 20;
    
    size_t numIteration = 1000;
    
    size_t m = 5;
    
    size_t dim = 3;
    
    std::string filename = "config";
    
    void setHistorySize(const size_t hs) { m = hs; }
    
    void setNumFile(const size_t nf) { numFile = nf; }
    
    void setNumIteration(const size_t i) { numIteration = i; }
    
    void setFileName(const std::string fn) { filename = fn; }
    
    void setDim(const size_t d) { dim = d; }
    
    
    void minimize(Problem<T> &objFunc, Vector<T> & x0) {
        
        std::ofstream outputfile;
    
        const size_t DIM = x0.rows();
        double minimum_energy = objFunc.value(x0);
        printf("\n\n\n");
        
        
        const size_t NUMPTS = DIM / dim;

        Matrix<T> sVector = Matrix<T>::Zero(DIM, m);
        Matrix<T> yVector = Matrix<T>::Zero(DIM, m);

        Vector<T> alpha = Vector<T>::Zero(m);
        Vector<T> grad(DIM), q(DIM), grad_old(DIM), s(DIM), y(DIM), minimum_x(DIM);
        objFunc.gradient(x0, grad);
        Vector<T> x_old = x0;

        size_t iter = 0;
        double H0k = 1;
        double alpha_init = 1;
        double md;
        double mg;
        int tmp = numIteration/numFile;
        

        this->m_current.reset();
        do {
            const T relativeEpsilon = static_cast<T>(0.0001) * std::max(static_cast<T>(1.0), x0.norm());

            if (grad.norm() < relativeEpsilon)
                break;

            //Algorithm 7.4 (L-BFGS two-loop recursion)
            q = grad;
            const int k = std::min(m, iter);
            
            // for i = k − 1, k − 2, . . . , k − m
            for (int i = k - 1; i >= 0; i--) {
                // alpha_i <- rho_i*s_i^T*q
                double product = static_cast<Vector<T>>(sVector.col(i))
                .dot(static_cast<Vector<T>>(yVector.col(i)));
                const double rho = 1.0 / product;
                alpha(i) = rho * static_cast<Vector<T>>(sVector.col(i)).dot(q);
                // q <- q - alpha_i*y_i
                q = q - alpha(i) * yVector.col(i);
            }
            // r <- H_k^0*q
            
            q = H0k * q;
            
            //for i k − m, k − m + 1, . . . , k − 1
            for (int i = 0; i < k; i++) {
                // beta <- rho_i * y_i^T * r
                double product = static_cast<Vector<T>>(sVector.col(i))
                .dot(static_cast<Vector<T>>(yVector.col(i)));
                const double rho = 1.0 / product;
                const double beta = rho * static_cast<Vector<T>>(yVector.col(i)).dot(q);
                // r <- r + s_i * ( alpha_i - beta)
                q = q + sVector.col(i) * (alpha(i) - beta);
            }
            
            // stop with result "H_k*f_f'=q"

            // any issues with the descent direction ?
            double descent = -grad.dot(q);
            
//            if (iter == 0)
                alpha_init =  1.0 / grad.norm();
            

            
            
            if (descent > -0.0001 * relativeEpsilon) {
                q = -1 * grad;
                iter = 0;
                alpha_init = 1.0;
            }
            
            printf("alpha_init: %f\n" , alpha_init);

            // find steplength
            Armijo<T, decltype(objFunc), 1>::linesearch(x0, -q,  objFunc, alpha_init) ;
            
            

            
        
            grad_old = grad;
            objFunc.gradient(x0, grad);

            s = x0 - x_old;
            y = grad - grad_old;

            // update the history
            if (iter < m) {
                sVector.col(iter) = s;
                yVector.col(iter) = y;
            } else {
                sVector.leftCols(m - 1) = sVector.rightCols(m - 1).eval();
                sVector.rightCols(1) = s;
                yVector.leftCols(m - 1) = yVector.rightCols(m - 1).eval();
                yVector.rightCols(1) = y;
            }
            
            
            
            
            
            
            double dot =  static_cast<double>(y.dot(y));
            if (dot <= 1e-7) {
                break;
//                double pp = objFunc.value(x0);
//                if (pp < minimum_energy) {
//                    minimum_energy = pp;
//                    minimum_x = x0;
//                }
//                printf("reach one local minimum, energy: %f\n\n\n", pp);
//                alpha_init = 1/grad.norm();
//                H0k = 1;
            }
            else{
                // update the scaling factor
                H0k = y.dot(s) / dot;
                
                // update alpha_init
                objFunc.findInitial(md, mg);
                alpha_init = 0.5 / (mg*NUMPTS);
            }
            
            
            x_old = x0;
            
            
            // std::cout << "iter: "<<globIter<< ", f = " <<  objFunc.value(x0) << ", ||g||_inf "
            // <<gradNorm  << std::endl;
            iter++;
            
            std::cout << "iter: " << iter << std::endl;
            
            
            if (iter % tmp == 0) {
                std::string name = "./output/" + filename + std::to_string(iter/tmp) + ".txt";
                outputfile.open(name.c_str());
                if (outputfile.fail()) {
                    printf("Error writing file\n");
                }
                
                writeFiles(outputfile, name, x0, 3);
                
                outputfile.close();
            }
            
            
            ++this->m_current.iterations;
            this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
            this->m_status = checkConvergence(this->m_stop, this->m_current);
        } while (this->m_status == Status::Continue && iter < numIteration);

        if (minimum_x.norm() > 1e-5) {
            x0 = minimum_x;
        }
        
        
    }

};

}
/* namespace cppoptlib */

#endif /* LBFGSSOLVER_H_ */
