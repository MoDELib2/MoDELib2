/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VTKsegments_cpp_
#define model_VTKsegments_cpp_


#include <VTKsegments.h>

namespace model
{

template <int dim>
double getEwaldLength(const std::vector<Eigen::Matrix<double,dim,1>>& periodicBasis,const double& EwaldLengthFactor)
{
    
    if(periodicBasis.size())
    {
        // Compute generalized volume of periodic lattice cell
        Eigen::MatrixXd B(Eigen::MatrixXd::Zero(3,periodicBasis.size()));
        for(size_t k=0;k<periodicBasis.size();++k)
        {
            B.col(k)=periodicBasis[k];
        }
        const double vol(sqrt((B.transpose()*B).determinant())/CTM::factorial(periodicBasis.size()));
//        std::cout<<"vol="<<vol<<std::endl;
//        std::cout<<"edge="<<std::pow(vol,1.0/periodicBasis.size())<<std::endl;
//        std::cout<<"elength="<<EwaldLengthFactor*std::pow(vol,1.0/periodicBasis.size())<<std::endl;
        return EwaldLengthFactor*std::pow(vol,1.0/periodicBasis.size());
    }
    else
    {
        return 0.0;
    }
}

    VTKsegments::VTKsegments(const std::string& folderName) :
    ///* init */ traitsIO(folderName)
    /* init */ perser(folderName+"/inputFiles/vtkSegments.txt")
    /* init */,material(perser.readString("materialFile",true),perser.readScalar<double>("absoluteTemperature",true))
    /* init */,quadPerLength(perser.readScalar<double>("quadPerLength",true))
    /* init */,externalStress(perser.readMatrix<double>("ExternalStress0",3,3,true))
    /* init */,stochasticForceGenerator(nullptr)
    /* init */,C2G(perser.readMatrix<double>("C2G",3,3,true))
    {
        
    }

//    typename VTKsegments::VectorDimI VTKsegments::getPbcFlags(const std::string& filename) const
//    {
//
//    }


    const std::vector<DislocationQuadraturePoint<3,0>>& VTKsegments::quadraturePoints() const
    {
        return *this;
    }

    std::vector<DislocationQuadraturePoint<3,0>>& VTKsegments::quadraturePoints()
    {
        return *this;
    }

    // const std::vector<StressStraight<3>>& VTKsegments::segments() const
    // {
    //     return *this;
    // }

    // std::vector<StressStraight<3>>& VTKsegments::segments()
    // {
    //     return *this;
    // }


    const std::map<std::pair<int,int>,StressStraight<3>> &VTKsegments::segments() const
    {
        return *this;
    }

    std::map<std::pair<int,int>,StressStraight<3>> &VTKsegments::segments()
    {
        return *this;
    }

    // const std::vector<typename VTKsegments::VectorDim>& VTKsegments::nodes() const
    // {
    //     return *this;
    // }

    const std::deque<typename VTKsegments::VectorDim> &VTKsegments::nodes() const
    {
        return *this;
    }


    std::deque<typename VTKsegments::VectorDim> &VTKsegments::nodes()
    {
        return *this;
    }

double VTKsegments::latticeParameter() const
{
    if(material.crystalStructure=="BCC")
    {
        return 2.0*material.b_SI/sqrt(3.0);
    }
    else if(material.crystalStructure=="FCC")
    {
        return 2.0*material.b_SI/sqrt(2.0);
    }
    else
    {
        std::cout<<"Unknown lattice parameter for "<<material.crystalStructure<<"'. Exiting."<<std::endl;
        exit(EXIT_FAILURE);
        return 0.0;
    }
}

// std::vector<typename VTKsegments::VectorDim>& VTKsegments::nodes()
// {
//     return *this;
// }


    void VTKsegments::writeVTK(const std::string& vtkFilePrefix,std::map<std::pair<int,int>,std::pair<size_t,std::set<int>>> &segIDmap) const
    {
        const std::string quadFileName(vtkFilePrefix+"_quadrature.vtk");
        std::ofstream quadFile(quadFileName);
        quadFile<<"# vtk DataFile Version 3.0\n";
        quadFile<<"# Dislocation lines converted from MoDELib file\n";
        quadFile<<"ASCII\n";
        quadFile<<"DATASET UNSTRUCTURED_GRID\n";
        quadFile<<"POINTS "+std::to_string(quadraturePoints().size()+nodes().size())+" double\n";

        for(const auto& node : nodes())
        {// write all nodes positions
            quadFile<<std::setprecision(15)<<std::scientific<<node.transpose()*material.b_SI*1.0e10<<"\n";
        }
        
        for(const auto& qp : quadraturePoints())
        {// write all nodes positions
            quadFile<<std::setprecision(15)<<std::scientific<<qp.r.transpose()*material.b_SI*1.0e10<<"\n";
        }
        
        // Create map of quadraturePoints by segments
        std::map<std::pair<size_t,size_t>,std::set<size_t>> qPointMap;
        for(size_t q=0;q<quadraturePoints().size();++q)
        {
            const auto& qp(quadraturePoints()[q]);
            const std::pair<size_t,size_t> key(std::make_pair(qp.sourceID,qp.sinkID));
            qPointMap[key].insert(q);
        }
        
        if(qPointMap.size()!=segments().size())
        {
            std::cout<<"qPointMap.size()="<<qPointMap.size()<<std::endl;
            std::cout<<"segments().size()="<<segments().size()<<std::endl;
            throw std::runtime_error("qPointMap.size()!=segments().size()");
        }
        
        quadFile<<"\nCELLS "+std::to_string(qPointMap.size())+" "+std::to_string(quadraturePoints().size()+3*qPointMap.size())+"\n";
        for(const auto& pair : qPointMap)
        {
            quadFile<<pair.second.size()+2<<" "<<pair.first.first<<" ";
            for(const auto& q : pair.second)
            {
                quadFile<<q+nodes().size()<<" ";
            }
            quadFile<<pair.first.second<<"\n";
        }
        
         quadFile << "\nCELL_TYPES " + std::to_string(segments().size()) + "\n";
        for (const auto &pair : qPointMap)
        {
            quadFile << 4 << "\n";
        }
        quadFile << "\nPOINT_DATA " + std::to_string(quadraturePoints().size() + nodes().size()) + "\n";
        quadFile << "VECTORS PK_force double\n";
        for (const auto &node : nodes())
        {
            quadFile << std::setprecision(15) << std::scientific << Eigen::Matrix<double, 1, 3>::Zero() << "\n";
        }
        for (const auto &outerPair : qPointMap)
        {
            for (const auto &innerPair : outerPair.second)
            {
                const auto &qp(quadraturePoints()[innerPair]);
                quadFile << std::setprecision(15) << std::scientific << qp.pkForce.transpose() * material.b_SI * material.mu_SI * 10.0 / 160.21766208 << "\n";
                //quadFile << std::setprecision(15) << std::scientific << qp.pkForce.transpose()  << "\n";
            }
        }
        quadFile << "\nTENSORS stress double\n";

        for (const auto &node : nodes())
        {
            quadFile << std::setprecision(15) << Eigen::Matrix<double, 3, 3>::Zero() << "\n\n";
        }

        for (const auto &outerPair : qPointMap)
        {
            for (const auto &innerPair : outerPair.second)
            {
                const auto &qp(quadraturePoints()[innerPair]);
                quadFile << std::setprecision(15) << std::scientific << qp.stress * material.mu_SI * 1e-9 << "\n\n";
            }
        }
        quadFile << "\nCELL_DATA " + std::to_string(qPointMap.size()) + "\n";
        quadFile << "SCALARS dislocation_index int\n";
        quadFile << "LOOKUP_TABLE default\n";
        for (size_t k = 0; k < qPointMap.size(); ++k)
        {
            quadFile << k << "\n"; // not sure if this is the right number to write
        }
        //NEED TO ENABLE THIS SECTION
        quadFile << "\nVECTORS burgers_vector_global double\n";

        for (const auto &pair : segIDmap)
        {
        
            const auto segID = pair.first;
            const auto &seg(segments().at(segID));
            std::cout<<"SegID: "<<segID.first<<"-> "<<segID.second <<"-> seg.b: "<<seg.b.transpose() * material.b_SI * 1.0e10<<std::endl;
            quadFile << seg.b.transpose() * material.b_SI * 1.0e10 << "\n";
        }

        quadFile << "\nVECTORS burgers_vector_local double\n";
        const auto G2C(C2G.inverse());
        std::cout<<"C2G: "<<C2G<<std::endl;
        std::cout<<"G2C: "<<G2C<<std::endl;

        for (const auto &pair : segIDmap)
        {
        
            const auto segID = pair.first;
            const auto &seg(segments().at(segID));
            
            quadFile << (G2C * seg.b).transpose() * material.b_SI / latticeParameter() << "\n";
        
        }


    }

    void VTKsegments::ExternalController()
    {
       const auto externalStress=perser.readMatrix<double>("externalStress",3,3,true);
    }

    void VTKsegments::readCAFile(const std::string& vtkFilePrefix)
    {
        VectorDimI pbcFlags;
        const std::string caFileName(vtkFilePrefix+".ca");
        std::ifstream caFile(caFileName); //access vtk file
        
        
        
        if(caFile.is_open())
        {
            
            std::string line;
            while (std::getline(caFile, line)) //begin parsing vtk file for lines
            {
                if(line.find("SIMULATION_CELL_MATRIX")!=std::string::npos) //if POINTS is read, read npos
                {
                    for(int d=0;d<dim;++d)
                    {
                        std::getline(caFile, line);
                        std::stringstream ss(line);//store lines of vtk file in ss
                        ss>>cellMatrix(d,0)>>cellMatrix(d,1)>>cellMatrix(d,2);//store nodal positions to x,y,z
                    }
                    std::cout<<"SIMULATION_CELL_MATRIX=\n"<<cellMatrix<<std::endl;
                    
                    std::vector<Eigen::Matrix<double,dim,1>> cellVectors;
                    for(int d=0;d<3;++d)
                    {
                        cellVectors.emplace_back(cellMatrix.row(d)*1.0e-10/material.b_SI);
                    }
                    EwaldLength=getEwaldLength(cellVectors,perser.readScalar<double>("EwaldLengthFactor",true));
                    std::cout<<"Ewald Length"<<std::endl;
                    std::cout<<EwaldLength<<std::endl;
                }
                
                if(line.find("PBC_FLAGS")!=std::string::npos) //if POINTS is read, read npos
                {
                    std::string temp;
                    std::stringstream ss(line);//store lines of vtk file in ss
                    ss>>temp;
                    ss>>pbcFlags(0)>>pbcFlags(1)>>pbcFlags(2);
                    std::cout<<"PBC_FLAGS="<<pbcFlags.transpose()<<std::endl;

                }
                
            }
            
            periodicShifts.clear();
            for(int i=-pbcFlags(0);i<pbcFlags(0)+1;++i)
            {
                for(int j=-pbcFlags(1);j<pbcFlags(1)+1;++j)
                {
                    for(int k=-pbcFlags(2);k<pbcFlags(2)+1;++k)
                    {
                        periodicShifts.push_back( i*cellMatrix.row(0)*1.0e-10/material.b_SI
                                                 +j*cellMatrix.row(1)*1.0e-10/material.b_SI
                                                 +k*cellMatrix.row(2)*1.0e-10/material.b_SI);
                    }
                }
            }
            std::cout<<"periodicShifts.size()="<<periodicShifts.size()<<std::endl;
            
            //prints PBCvectors
            std::cout<<"PBCvectors"<<std::endl;
            for (const auto& matrix : periodicShifts) 
            {
            std::cout << matrix << std::endl;
            }
            
        }
        else
        {
            throw std::runtime_error("Cannot open file "+ caFileName);
        }
        
    }


    void VTKsegments::updateQuadraturePoints(const std::string& vtkFilePrefix, MatrixDim externalStress, std::map<std::pair<int,int>,std::pair<size_t,std::set<int>>> &segIDmap )
    {
        
        
        
        
        
        
    #ifdef _OPENMP
        const size_t nThreads = omp_get_max_threads();
    #else
        const size_t nThreads = 1;
    #endif
        std::cout <<"updating QuadraturePoints (" << nThreads << " threads) " << std::flush;
        const auto t0= std::chrono::system_clock::now();

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
        for(size_t q=0;q<quadraturePoints().size();++q)
        {
            auto& qPoint(quadraturePoints()[q]);
            for(const auto& seg : segments())
            {
                for(const auto& shift : periodicShifts)
                {
                    qPoint.stress+=seg.second.stress(qPoint.r+shift);
                }
            }
        }
        for (size_t q = 0; q < quadraturePoints().size(); ++q)
        {
            auto &qPoint(quadraturePoints()[q]);
            qPoint.stress += externalStress;
        }

        for (const auto &pair : segIDmap)
        {
            const auto segID = pair.first;
            const auto &seg(segments().at(segID));
            for (const auto &qID : pair.second.second)
            {
                auto &qPoint(quadraturePoints()[qID]);
                qPoint.pkForce = (qPoint.stress * seg.b).cross(qPoint.rl);
            }
        }

        std::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

        writeVTK(vtkFilePrefix,segIDmap);
    }

    void VTKsegments::readVTK(const std::string& vtkFilePrefix)
    {
        nodes().clear();
        segments().clear();
        quadraturePoints().clear();
        
        const std::string vtkFileName(vtkFilePrefix+".vtk");
        std::ifstream vtkFile(vtkFileName); //access vtk file
        std::map<std::pair<int,int>, std::pair<size_t,std::set<int>>> segIDmap; // key= [segID] value=[cellID, set of qID's]
        

        readCAFile(vtkFilePrefix);

        if(vtkFile.is_open())
        {
//            std::vector<VectorDim> rawPoints;
            std::vector<std::vector<size_t>> cells;
            std::vector<VectorDim> rawBurgers;
            const auto sfCoeffs(SplineSegmentBase<3,0>::sfCoeffs(0.0)); // shape function coefficients for linear segments
            
            int loopsSize(0);
            std::vector<VectorDim> nAVector;
            std::string line;
            
            
            while (std::getline(vtkFile, line)) //begin parsing vtk file for lines
            {
                if(line.find("POINTS")!=std::string::npos) //if POINTS is read, read npos
                {
                    const size_t firstSpace(line.find(' ')); //vtk file formatting const
                    const size_t secondSpace(line.find(' ',firstSpace+1)); //vtk file formatting const
                    
                    if(secondSpace>firstSpace) // if the line indentation is longer...
                    {
                        const size_t nodesSize(std::atoi(line.substr(firstSpace+1,secondSpace-firstSpace-1).c_str())); //read the number of nodes
                        
                        std::cout<<"Reading "<<nodesSize<<" nodes"<<std::endl; //print the number of nodes
                        
                        double x,y,z; //create doubles to store node positions
                        
                        for(size_t n=0;n<nodesSize;++n) //index over all nodes
                        {
                            std::getline(vtkFile, line); //access lines of vtk file
                            std::stringstream ss(line);//store lines of vtk file in ss
                            ss>>x>>y>>z;//store nodal positions to x,y,z
                            nodes().push_back((VectorDim()<<x,y,z).finished()*1.0e-10/material.b_SI);//construct nodes()  [<x_1,x_2,x_3>,loopNo]
                        }
                    }
                    else
                    {//POINTS was not read
                        throw std::runtime_error("Unable to extract number of nodes from line: "+line);
                    }
                }
                
                if(line.find("CELLS")!=std::string::npos) //if CELLS is read, read npos
                {
                    const size_t firstSpace(line.find(' ')); //vtk file formatting const
                    const size_t secondSpace(line.find(' ',firstSpace+1)); //vtk file formatting const
                    
                    if(secondSpace>firstSpace) // if the line indentation is longer...
                    {
                        
                        loopsSize=std::atoi(line.substr(firstSpace+1,secondSpace-firstSpace-1).c_str()); //read the total number of loops loopsSize
                        
                        std::cout<<"Reading "<<loopsSize<<" cells"<<std::endl; // print the total number of loops
                        //                    cells.emplace_back();
                        for(int n=0;n<loopsSize;++n) // index over the number of loops
                        {
                            std::getline(vtkFile, line); // read the contents of cells
                            std::stringstream ss(line); // store the contents of cells in ss
                            cells.emplace_back(); //populate cells
                            int temp; //creante int temp
                            int count =0; //create a count initialized at 0
                            
                            while (ss >> temp) // while the nodal index is less than temp...
                            {
                                if(count==0) //if count = 0
                                {
                                    count++; // increase count by 1
                                    continue; //rerun the loop
                                }
                                cells.back().push_back(temp); //populate cells wih node connections <node pts, Burgers vector>
                            }
                        }
                    }
                    else
                    {//CELLS was not read
                        throw std::runtime_error("Unable to extract number of loops from line: "+line);
                    }
                }
                
                
                if(line.find("burgers_vector_world")!=std::string::npos) //if burgers_vector_world is read, read npos
                {
                    double x,y,z; //initialize x,y,z, for burgers vector components
                    
                    for(int n=0;n<loopsSize;++n) //index over the number of loops
                    {
                        std::getline(vtkFile, line);//read the contents of burgers_vector_world
                        std::stringstream ss(line);// store the contents of burgers_vector_world in ss
                        ss>>x>>y>>z;// pass ss to x,y,z
                        rawBurgers.push_back((VectorDim()<<x,y,z).finished()*1.0e-10/material.b_SI);
                    }
                }
            }
            
            if(rawBurgers.size()==cells.size())
            {

                for (size_t k = 0; k < cells.size(); ++k)
                {
                    const auto Burgers(rawBurgers[k]);
                    VectorDim nA(VectorDim::Zero()); // right-handed loop normal
                    const VectorDim P0(nodes()[0]); //sample point P0 from vtk file
                   

                    for(int n=0;n<cells[k].size();++n)//index over all nodes
                    {
                        const int n1(n<cells[k].size()-1? n+1 : 0);//index k1 to be one less than the current node number

                        //std::cout<<"k1= "<<k1<<std::endl;// print the value of k1

                        const size_t sourceID(cells[k][n]);
                        const size_t sinkID(cells[k][n1]);
                        const VectorDim &sourceP(nodes()[sourceID]);
                        const VectorDim &sinkP(nodes()[sinkID]);
                        //std::cout<<"sourcePos->sinkPos: "<<sourcePos<<" -> "<<sinkPos<<std::endl;

                        const double linkNorm((sinkP-sourceP).norm());//returns the normal to the connection

                        
                        nA+=0.5*(sourceP-P0).cross(sinkP-sourceP);//modify the right handed loop normal
                    }

                    nAVector.push_back(nA.normalized()); //normalize the right handed loop normal

                    for (size_t n = 0; n < cells[k].size() - 1; ++n)
                    {
                        const size_t sourceID(cells[k][n]);
                        const size_t sinkID(cells[k][n + 1]);

                        
                        if (sinkID>sourceID)
                        {
                            const VectorDim &sourceP(nodes()[sourceID]);
                            const VectorDim &sinkP(nodes()[sinkID]);
                            const std::pair<int,int> sourceSink=std::make_pair(sourceID,sinkID);
                            
                            // std::cout<<"Ewald1"<<std::endl;
                            // std::cout<<EwaldLength<<std::endl;

                            const std::pair<std::pair<int,int>,StressStraight<3>> temp=std::make_pair(sourceSink,StressStraight<3>(material, sourceP, sinkP, Burgers, EwaldLength));
                            segments().insert(temp
                                        //std::piecewise_construct,
                                        //std::make_tuple(sourceID,sinkID),
                                        //std::make_tuple(material, sourceP, sinkP, Burgers)
                                        );


                            // segments().emplace(std::piecewise_construct,
                            //             std::make_tuple(sourceID,sinkID),
                            //             std::make_tuple(material, sourceP, sinkP, Burgers, EwaldLength)
                            //             );

                            // const VectorDim chord(sinkP - sourceP);
                            // const double chordLength(chord.norm());
                            // const int qOrder = QuadratureDynamicType::lowerOrder(quadPerLength * chordLength);
                            // const MatrixNcoeffDim dofs((MatrixNcoeffDim() << sourceP.transpose(), sinkP.transpose()).finished());

                            // for (int q = 0; q < qOrder; ++q)
                            // {

                            //     quadraturePoints().emplace_back(sourceID, sinkID, q, qOrder, sfCoeffs, dofs);
                                
                                segIDmap[std::make_pair(sourceID,sinkID)].first=k;
                            //     segIDmap[std::make_pair(sourceID,sinkID)].second.insert(quadraturePoints().size() - 1);
                            //     //   segIDmap[std::pair(k,std::make_pair(sourceID,sinkID))].insert(quadraturePoints().size() - 1);
                            // }

                        }
                        else
                        {
                            const VectorDim &sourceP(nodes()[sinkID]);
                            const VectorDim &sinkP(nodes()[sourceID]);
                            const std::pair<int,int> sourceSink=std::make_pair(sinkID,sourceID);
                            // std::cout<<"Ewald2"<<std::endl;
                            // std::cout<<EwaldLength<<std::endl;

                            const std::pair<std::pair<int,int>,StressStraight<3>> temp=std::make_pair(sourceSink,StressStraight<3>(material, sourceP, sinkP, -Burgers, EwaldLength));
                                segments().insert(temp
                                        //std::piecewise_construct,
                                        //std::make_tuple(sinkID,sourceID),
                                        //std::make_tuple(material, sourceP, sinkP, -Burgers)
                                        );



                            //  segments().emplace(std::piecewise_construct,
                            //             std::make_tuple(sinkID,sourceID),
                            //             std::make_tuple(material, sourceP, sinkP, -Burgers, EwaldLength)
                            //             );

                             segIDmap[std::make_pair(sinkID,sourceID)].first=k;      
                           
                        }
                        //const std::pair<int,int> key=std::make_pair(std::min(sourceID,sinkID),std::max(sourceID,sinkID));
                        

                        //segments().emplace_back(material, sourceP, sinkP, Burgers);
                        
                    }
                }
                for (const auto& seg : segments())
                {
                    const double chordLength(seg.second.length);
                    const int sourceID=seg.first.first;
                    const int sinkID=seg.first.second;
                    const VectorDim sourceP=seg.second.P0;
                    const VectorDim sinkP=seg.second.P1;
                    const int qOrder = QuadratureDynamicType::lowerOrder(quadPerLength * chordLength);
                    const MatrixNcoeffDim dofs((MatrixNcoeffDim() << sourceP.transpose(), sinkP.transpose()).finished());

                    for (int q = 0; q < qOrder; ++q)
                    {

                        quadraturePoints().emplace_back(sourceID, sinkID, q, qOrder, sfCoeffs, dofs);
                        
                        //segIDmap[std::make_pair(sourceID,sinkID)].first=k;
                        segIDmap[std::make_pair(sourceID,sinkID)].second.insert(quadraturePoints().size() - 1);
                        //segIDmap[std::pair(k,std::make_pair(sinkID,sourceID))].insert(quadraturePoints().size() - 1);
                    }
                }
              
                std::cout << "segments().size()=" << segments().size() << std::endl;
                std::cout << "quadraturePoints().size()=" << quadraturePoints().size() << std::endl;
            }
            else
            {
                std::cout << "rawBurgers.size()=" << rawBurgers.size() << std::endl;
                std::cout << "cells.size()=" << cells.size() << std::endl;
                throw std::runtime_error("rawBurgers.size() NOT EQUAL TO cells.size()");
            }
            
        }
        else
        {
            std::cout<<"Cannot open "<<vtkFileName<<std::endl;
        }
        
        updateQuadraturePoints(vtkFilePrefix, externalStress, segIDmap);
        
    }

}
#endif
