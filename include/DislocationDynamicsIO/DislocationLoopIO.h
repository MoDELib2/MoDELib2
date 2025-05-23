/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopIO_H_
#define model_DislocationLoopIO_H_

#include <iomanip>
#include <tuple>
#include <Eigen/Dense>

namespace model
{
    
    template<short unsigned int dim>
    struct DislocationLoopIO
    {
        
        enum DislocationLoopType{GLISSILELOOP,SESSILELOOP,VIRTUALLOOP};

        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
         size_t sID;          // sID
         VectorDim B;          // position
         VectorDim N;          // normal
         VectorDim P;          // position
         size_t grainID;          // component ID
         int loopType;
         std::tuple<double,double,double,double> loopLength;
         double slippedArea;
        
        template<typename DislocationLoopType>
        DislocationLoopIO(const DislocationLoopType& dL) :
        /* init */ sID(dL.sID)
        /* init */,B(dL.flow().cartesian())
        /* init */,N(dL.glidePlane? (dL.slippedArea()>FLT_EPSILON? dL.rightHandedUnitNormal() : dL.glidePlane->unitNormal) : VectorDim::Zero())
        /* init */,P(dL.glidePlane? dL.glidePlane->P : (*dL.loopLinks().begin())->source->get_P() )
        /* init */,grainID(dL.grain.grainID)
        /* init */,loopType(dL.loopType)
        /* init */,loopLength(dL.network().outputLoopLength? dL.loopLength() : std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(dL.slippedArea())
        {
            
        }
        
        DislocationLoopIO(const size_t& sID_in,         // sID
                          const VectorDim& B_in,          // position
                          const VectorDim& N_in,          // velocity
                          const VectorDim& P_in,          // velocity
                          const size_t& grainID_in,
                          const int& loopType_in
                          ) :
        /* init */ sID(sID_in)
        /* init */,B(B_in)
        /* init */,N(N_in)
        /* init */,P(P_in)
        /* init */,grainID(grainID_in)
        /* init */,loopType(loopType_in)
        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(0.0)
        {
            
        }
        
        DislocationLoopIO() :
        /* init */ sID(0)
        /* init */,B(VectorDim::Zero())
        /* init */,N(VectorDim::Zero())
        /* init */,P(VectorDim::Zero())
        /* init */,grainID(0)
        /* init */,loopType(0)
        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(0.0)
        {
            
        }
        
        DislocationLoopIO(std::stringstream& ss) :
        /* init */ sID(0)
        /* init */,B(VectorDim::Zero())
        /* init */,N(VectorDim::Zero())
        /* init */,P(VectorDim::Zero())
        /* init */,grainID(0)
        /* init */,loopType(0)
        /* init */,loopLength(std::make_tuple(0.0,0.0,0.0,0.0))
        /* init */,slippedArea(0.0)
        {
            ss>>sID;
            for(int d=0;d<dim;++d)
            {
                ss>>B(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>N(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>P(d);
            }
            ss>>grainID;
            ss>>loopType;
            double l1,l2,l3,A;
            ss>>l1>>l2>>l3>>A;
            loopLength=std::make_tuple(l1,l2,l3,A);
            ss>>slippedArea;
        }
        
        template <class T>
        friend T& operator << (T& os, const DislocationLoopIO<dim>& ds)
        {
            os  << ds.sID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.B.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.N.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ds.P.transpose()<<"\t"
            /**/<< ds.grainID<<"\t"
            /**/<< ds.loopType<<"\t"
            /**/<< std::get<0>(ds.loopLength)<<"\t"<< std::get<1>(ds.loopLength)<<"\t"<< std::get<2>(ds.loopLength)<<"\t"<< std::get<3>(ds.loopLength)<<"\t"
            /**/<< ds.slippedArea;
            return os;
        }
        
	};
	
}
#endif

