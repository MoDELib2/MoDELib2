/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_GammaSurface_cpp_
#define model_GammaSurface_cpp_

#include <memory>
#include <assert.h>
#include <TerminalColors.h>
#include <LatticeModule.h>
#include <GammaSurface.h>

namespace model
{

//    typename GammaSurface::MatrixDim GammaSurface::getG2L(const VectorDim& x,
//                                                          const VectorDim& z)
//    {
//        const double xNorm(x.norm());
//        const double zNorm(z.norm());
//        assert(xNorm>FLT_EPSILON);
//        assert(zNorm>FLT_EPSILON);
//        assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
//        Eigen::Matrix3d temp(Eigen::Matrix3d::Identity());
//        temp.col(2)=z/zNorm;
//        temp.col(0)=x/xNorm;
//        temp.col(1)=temp.col(2).cross(temp.col(0));
//        return temp.transpose();
//    }
//
//GammaSurface::MatrixLowerDim GammaSurface::getLocalBasis(const LatticeVector<dim>& a1,
//                                                         const LatticeVector<dim>& a2)
//{
//    const Eigen::Matrix3d R(getG2L(a1.cartesian(),a1.cross(a2).cartesian()));
//    MatrixLowerDim temp(MatrixLowerDim::Zero());
//    temp.col(0)=(R*a1.cartesian()).segment<2>(0);
//    temp.col(1)=(R*a2.cartesian()).segment<2>(0);
//    return temp;
//}


GammaSurface::GammaSurface(const MatrixLowerDim& A,
                           const Eigen::Matrix<double,Eigen::Dynamic,lowerDim>& waveVectors,
                           const Eigen::Matrix<double,Eigen::Dynamic,lowerDim+1>& f,
                           const int& rotSymm,
                           const std::vector<Eigen::Matrix<double,lowerDim,1>>& mirSymm) :
/* init */ PeriodicLatticeInterpolant<2>(A,waveVectors,f,rotSymm,mirSymm)
///* init */,G2L(getG2L(a1.cartesian(),a1.cross(a2).cartesian()))
{
    std::cout<<greenColor<<"Creating GammaSurface with local basis\n "<<defaultColor<<this->A<<std::endl;
}

//    double GammaSurface::operator()(const VectorDim& b)
//    {
//        const VectorDim bL(G2L*b);
//        if(std::fabs(bL(dim-1))>FLT_EPSILON)
//        {
//            throw std::runtime_error("SLIP VECTOR NOT ON GAMMA-SURFACE PLANE");
//        }
//        return PeriodicLatticeInterpolant<2>::operator()(bL.segment<lowerDim>(0));
//    }

}
#endif
