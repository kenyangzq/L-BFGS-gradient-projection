
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include "meta.h"
#include "problem.h"
#include "solver/lbfgssolver.h"
#include "solver/bfgssolver.h"
#include "solver/lbfgsbsolver.h"

using namespace std;



double dist(const cppoptlib::Vector<double> & p1,const cppoptlib::Vector<double> & p2){
    int dim = p1.rows();
    double sum = 0;
    for (int i = 0; i < dim; ++i) {
        sum += (p1(i)-p2(i))*(p1(i)-p2(i));
    }
    return sqrt(sum);
}


// cutoff functions
double cutoff(double distance, double cutoff_radius){
    return pow(1-pow(distance/cutoff_radius, 4), 3);
}

double cutoffPrime(double distance, double cutoff_radius){
    double t = distance / cutoff_radius;
    return 12*pow(1-pow(t,4),2)*pow(t,3)/distance/cutoff_radius;
}



class minimize : public cppoptlib::Problem<double> {
    
    double cutoff_radius;
    double s;
    int cubes_per_side;
    int dim;
    int numpts;
    cppoptlib::Matrix<double> Cubes;
    cppoptlib::Matrix<double> pts3D;
    double minimal_distance;
    double minimal_gradient;
    
    
public:
    
    minimize(double r, double s_value, int d, int n, int c, int c_cube, int max_neighbor)
    :cutoff_radius(r), s(s_value), cubes_per_side(c), dim(d), numpts(n), Cubes(max_neighbor+1, c_cube),
    pts3D(dim, numpts), minimal_distance(0), minimal_gradient(0){};
    
    
    
    double value(const cppoptlib::Vector<double> &x){

        double total_energy = 0;
        
        //int max_neighbor = Cubes.rows();
        cppoptlib::Vector<int> neighbor_cube_indices (pow(3, dim));

        // get the 3d points matrix -- pts
        // assign points to cubes
        // the first field in each column is the number of points inside that column.
        BuildIndex(x);
        
        // find the neighbor_cube_indices for each cube and calculate the energy.
        for (int i = 0; i < Cubes.cols(); ++i)
        {
            neighbor_cube_indices = cppoptlib::Vector<int>::Ones(pow(3, dim))*(-1);
            findNeighborCubes(i, neighbor_cube_indices);
            //cout << neighbor_cube_indices.transpose() << endl << endl;
            int points_in_cube_number = Cubes(0, i);
            // j goes over all points in the i-th cube
            
            for (int j = 1; j <= points_in_cube_number; ++j) {
                int point_index = Cubes(j, i); // absolute point index
                total_energy += truncatedEnergy(x, neighbor_cube_indices, point_index);
            }
        }
        cout << fixed;
        cout << "Energy: " << total_energy << endl;
        
        return total_energy;
    }
    
    
    
    
    
    void gradient(const cppoptlib::Vector<double> &x, cppoptlib::Vector<double> &grad) {
    
        int max_neighbor = Cubes.rows();
        double distance;
        cppoptlib::Vector<int> neighbor_cube_indices (pow(3, dim));
        cppoptlib::Vector<double> temp_sum(dim), temp(dim);
        temp_sum.setZero();
        BuildIndex(x);
        double mgradient = 1000000;
        
        
        for (int index_cube = 0; index_cube < Cubes.cols(); ++index_cube)
        {
            neighbor_cube_indices = cppoptlib::Vector<int>::Ones(pow(3, dim))*(-1);
            findNeighborCubes(index_cube, neighbor_cube_indices);
            
            int points_in_cube = Cubes(0, index_cube);
            for (int j = 1; j <= points_in_cube; ++j)
            {
                temp_sum.setZero();
                int point_index = Cubes(j, index_cube);
                for (int k = 0; k < neighbor_cube_indices.size(); ++k)
                {
                    int tmp = neighbor_cube_indices(k);
                    if (tmp != -1)
                    {
                        int points_in_other_cube = Cubes(0, tmp);
                        for (int l = 1; l <= points_in_other_cube; l++)
                        {
                            int other_point_index = Cubes(l, tmp);
                            temp = pts3D.col(point_index) - pts3D.col(other_point_index);
                            distance = sqrt(temp.dot(temp));
                            if (other_point_index != point_index && distance < cutoff_radius)
                            {
                                temp_sum += (cutoffPrime(distance, cutoff_radius)*pow(distance, -s) +
                                             (-s)*cutoff(distance,cutoff_radius)*pow(distance, -s-2))*temp;
                            }
                        }
                    }
                }
                
                grad.segment(point_index*dim, dim) = temp_sum.transpose();
                double tmp = temp_sum.norm();
                if (tmp < mgradient)
                    mgradient = tmp;
            }
        }

        minimal_gradient = mgradient;
    }
    
    
    
    void pullBack(cppoptlib::Vector<double> &x){
        cppoptlib::Vector<double> temp_vector(dim);
        
        for (int i = 0; i < numpts; ++i) {
            temp_vector = x.segment(dim*i, dim);
            double norm = temp_vector.norm();
            if (norm < 1e-7) {
                cout << temp_vector.transpose() << endl;
                cout << "wtf" << endl;
            }
            x.segment(dim*i, dim) = temp_vector / norm;
        }
    }
    
    
    
    
////// all the helping method in the class start here ///////
    
    
    
    void findInitial(double & mdistance, double & mgradient){
        mdistance = minimal_distance;
        mgradient = minimal_gradient;
    }
    
    // x now is a 3d points vectors
    void BuildIndex(const cppoptlib::Vector<double> & x){
        
        Cubes.row(0).setZero();
        cppoptlib::Vector<double> temp_vector(dim);
        
        for (int i = 0; i < numpts; ++i) {
            temp_vector = x.segment(dim*i, dim);
            pts3D.col(i) = temp_vector.transpose();
            int max_neighbor = Cubes.rows();
            
            
            int cube_index = PointToCube(pts3D.col(i));
            if (cube_index == 70) {
                cout << pts3D.col(i);
                printf("cube index: %d\n", cube_index);
            }
            
            
            int tmp_numpts = Cubes(0, cube_index) + 1;
            if (tmp_numpts >= max_neighbor) {
                cout << "Warning: exceeding maximum neighbor; ignore point " << i << endl;
            }else{
                Cubes(tmp_numpts, cube_index) = i;
                Cubes(0, cube_index) = tmp_numpts;
            }
        }
        
        
    }
    
    
    int PointToCube(const cppoptlib::Vector<double> & thisPt){
        int tmp = 0, output = 0;
        for(int i = 0; i < dim; i++){
            tmp = floor((thisPt(i) +1.0) / cutoff_radius)+1;
            output += tmp * pow(cubes_per_side, i);
        }
        return output;
    }
    
    
    
    void findNeighborCubes(int cube_index, cppoptlib::Vector<int> & neighbors){
        
        cppoptlib::Vector<int> current_index(dim);
        //neighbors = Eigen::MatrixXd::Constant(neighbors.size(), 1, -1);
        int tmp = cube_index;
        
        for (int i = 0; i < dim; ++i){
            current_index(i) = tmp % cubes_per_side;
            tmp = tmp / cubes_per_side;
        }
        
        for (int i = 0; i < pow(3, dim); i++){
            int step = i;
            int flag = 0, index = 0;
            // index is the index of neighboring cube in the matrix of cubes;
            
            for (int j = 0; j < dim && flag == 0; j++){
                int k = step%3 - 1 + current_index[j];
                step = step/3;
                if (k < 0 || k >= cubes_per_side) flag=1;
                index += k * pow(cubes_per_side, j);
            }
            if (flag == 0) neighbors[i] = index;
        }
    }
    
    
    double truncatedEnergy(const cppoptlib::Vector<double> & V,
                        const cppoptlib::Vector<int> & neighbors, const int index){
        double energy = 0;
        double mdistance = 2;
        
        
        for (int i = 0; i < neighbors.size(); ++i)
        {
            int tmp = neighbors(i); // goes over all neighbor cubes
            if (tmp != -1) // check if there is a neighbor in this direction
            {
                int points_in_cube_number = Cubes(0, tmp);
                for (int j = 1; j <= points_in_cube_number; j++) {
                    int point_index = Cubes(j, tmp);
                    double distance = dist(V.segment(index*dim, dim), V.segment(point_index*dim, dim));
                    if (point_index != index && distance < cutoff_radius)
                        energy += pow(distance, -s)*cutoff(distance, cutoff_radius);
                    if (distance < mdistance) {
                        mdistance = distance;
                    }
                }
            }
        }
        minimal_distance = mdistance;
        
        return energy;
    }

    
    
    
};
























