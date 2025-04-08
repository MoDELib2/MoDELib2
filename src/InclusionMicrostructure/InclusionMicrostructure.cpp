/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionMicrostructure_cpp_
#define model_InclusionMicrostructure_cpp_

#include <InclusionMicrostructure.h>

namespace model
{
    /*!
    * @brief Constructor for InclusionMicrostructure.
    * @param mc Reference to the MicrostructureContainer.
    */

    template<int dim>
    InclusionMicrostructure<dim>::InclusionMicrostructure(MicrostructureContainerType& mc) :
    /* init */ MicrostructureBase<dim>("InclusionMicrostructure",mc)
    /* init */,InclusionMicrostructureBase<dim>(this->microstructures.ddBase)
    {
    }

    /*!
    * @brief Initializes inclusion configuration from input data.
    * 
    * @param configIO Input data structure containing configuration.
    * @param f_file Output file stream for field data.
    * @param F_labels Output file stream for field labels.
    * 
    * @throws std::runtime_error if phaseID is invalid, nodes are outside mesh,
    *         inclusion spans multiple grains, or node count is insufficient.
    */

    template<int dim>
    void InclusionMicrostructure<dim>::initializeConfiguration(const DDconfigIO<dim>& configIO,const std::ofstream& f_file,const std::ofstream& F_labels)
    {
        this->lastUpdateTime=this->ddBase.simulationParameters.totalTime;
        
        this->eshelbyInclusions().clear();
        this->polyhedronInclusionNodes().clear();
        for(const auto& inclusion : configIO.sphericalInclusions())
        {
            //        std::cout<<"Creating spherical inclusion "<<inclusion.inclusionID<<std::endl;
            const std::pair<bool,const Simplex<dim,dim>*> searchPair(this->ddBase.mesh.search(inclusion.C));
            if(searchPair.first)
            {
                
                const auto& grain(this->ddBase.poly.grain(searchPair.second->region->regionID));
                if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                {
                    const auto secondPhase(grain.singleCrystal->secondPhases().at(inclusion.phaseID));
                    EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                    std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new SphericalInclusion<dim>(inclusion.C,inclusion.a,inclusion.eT,this->ddBase.poly.nu,this->ddBase.poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                    this->eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                }
                else
                {
                    throw std::runtime_error("phaseID does not exist in grain.");
                }
            }
        }
        
        for(const auto& piNode : configIO.polyhedronInclusionNodes())
        {
            this->polyhedronInclusionNodes().emplace(piNode.nodeID,piNode);
        }
        
        std::map<size_t,std::map<size_t,std::vector<size_t>>> faceMap;
        for(const auto& edge : configIO.polyhedronInclusionEdges())
        {
            const size_t& iID(edge.inclusionID);
            const size_t& fID(edge.faceID);
            const size_t& sourceID(edge.sourceID);
            faceMap[iID][fID].push_back(sourceID);
        }
        
        
        for(const auto& inclusion : configIO.polyhedronInclusions())
        {
            //        std::cout<<"Creating polyhedron inclusion "<<inclusion.inclusionID<<std::endl;
            
            const auto faceIter(faceMap.find(inclusion.inclusionID));
            if(faceIter!=faceMap.end())
            {
                const auto& faces(faceIter->second);
                //            std::cout<<"    #faces= "<<faces.size()<<std::endl;
                std::set<const PolyhedronInclusionNodeIO<dim>*> uniquePolyNodes;
                for(const auto& pair : faces)
                {
                    for(const auto& nodeID : pair.second)
                    {
                        uniquePolyNodes.emplace(&this->polyhedronInclusionNodes().at(nodeID));
                    }
                }
                //            std::cout<<"    #nodes= "<<uniquePolyNodes.size()<<std::endl;
                if(uniquePolyNodes.size()>=dim+1)
                {
                    // Find grain
                    std::set<size_t> grainIDs;
                    for(const auto& nodePtr : uniquePolyNodes)
                    {
                        const std::pair<bool,const Simplex<dim,dim>*> searchPair(this->ddBase.mesh.search(nodePtr->P));
                        if(searchPair.first)
                        {
                            grainIDs.insert(searchPair.second->region->regionID);
                        }
                        else
                        {
                            throw std::runtime_error("inclusion node outside mesh");
                        }
                    }
                    
                    // Add inclusion
                    if(grainIDs.size()==1)
                    {
                        const auto& grain(this->ddBase.poly.grain(*grainIDs.begin()));
                        if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                        {
                            const auto secondPhase(grain.singleCrystal->secondPhases().at(inclusion.phaseID));
                            EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                            std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new PolyhedronInclusion<dim>(this->polyhedronInclusionNodes(),faces,inclusion.eT,this->ddBase.poly.nu,this->ddBase.poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                            this->eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                        }
                        else
                        {
                            throw std::runtime_error("phaseID does not exist in grain.");
                        }
                    }
                    else
                    {
                        throw std::runtime_error("inclusion across grain boundaries");
                    }
                }
                else
                {
                    throw std::runtime_error("inclusion does not have enough nodes");
                }
            }
            else
            {
                throw std::runtime_error("inclusionID not found in faceMap");
            }
        }
        
    }

    template<int dim>
    void InclusionMicrostructure<dim>::solve()
    {

    }

    /**
    * @brief Updates the internal state of the inclusion microstructure.
    */

    template<int dim>
    void InclusionMicrostructure<dim>::updateConfiguration()
    {
        this->lastUpdateTime=this->ddBase.simulationParameters.totalTime;
        
    }

    /**
    * @brief Returns the current time step size.
    * @return The current time step size.
    */

    template<int dim>
    double InclusionMicrostructure<dim>::getDt() const
    {
        return this->ddBase.simulationParameters.dtMax;
    }

    /**
    * @brief Outputs the current microstructure to configuration and auxiliary I/O.
    * @param configIO Output structure for the microstructure configuration.
    * @param auxIO Output structure for auxiliary data.
    * @param f_file Output stream for field data.
    * @param F_labels Output stream for label data.
    */

    template<int dim>
    void InclusionMicrostructure<dim>::output(DDconfigIO<dim>& configIO,DDauxIO<dim>& auxIO,std::ofstream& f_file,std::ofstream& F_labels) const
    {
        for(const auto& node : this->polyhedronInclusionNodes())
        {
            configIO.polyhedronInclusionNodes().emplace_back(node.second);
        }
        
        // Store Eshelby Inclusions
        for(const auto& ei : this->eshelbyInclusions())
        {
            
            auto* sphericalDerived = dynamic_cast<SphericalInclusion<dim>*>(ei.second.get());
            if (sphericalDerived)
            {
                configIO.sphericalInclusions().emplace_back(*sphericalDerived);
            }
            
            auto* polyhedronDerived = dynamic_cast<PolyhedronInclusion<dim>*>(ei.second.get());
            if (polyhedronDerived)
            {
                configIO.polyhedronInclusions().emplace_back(*polyhedronDerived);
                for(const auto& face : polyhedronDerived->faces)
                {
                    for(size_t k=0;k<face.second.size();++k)
                    {
                        const size_t k1(k<face.second.size()-1? k+1 : 0);
                        configIO.polyhedronInclusionEdges().emplace_back(polyhedronDerived->sID,face.first,face.second[k].first,face.second[k1].first);
                    }
                }
            }
        }
    }

    /**
    * @brief Returns the inelastic displacement rate of the inclusion.
    * @param x The position vector.
    * @param node Pointer to the node.
    * @param ele Pointer to the element.
    * @param guess Pointer to the simplex.
    * @return The inelastic displacement rate of the inclusion (zero vector).
    */
    template<int dim>
    typename InclusionMicrostructure<dim>::VectorDim InclusionMicrostructure<dim>::inelasticDisplacementRate(const VectorDim& x, const NodeType* const node, const ElementType* const ele,const SimplexDim* const guess) const
    {
        return VectorDim::Zero();
    }

    /**
    * @brief Returns the average plastic distortion of the inclusion.
    * @return The average plastic distortion of the inclusion (zero matrix).
    */

    template<int dim>
    typename InclusionMicrostructure<dim>::MatrixDim InclusionMicrostructure<dim>::averagePlasticDistortion() const
    {
        return MatrixDim::Zero();
    }

    /**
    * @brief Returns the average plastic distortion rate of the inclusion.
    * @return The average plastic distortion rate of the inclusion (zero matrix).
    */

    template<int dim>
    typename InclusionMicrostructure<dim>::MatrixDim InclusionMicrostructure<dim>::averagePlasticDistortionRate() const
    {
        return MatrixDim::Zero();
    }

    /**
    * @brief Returns the displacement of the inclusion.
    * @param VectorDimReference The position vector.
    * @param nodeType Type of Node.
    * @param ElementType Type of Element.
    * @param simplexType Type of Simplex.
    */

    template<int dim>
    typename InclusionMicrostructure<dim>::VectorDim InclusionMicrostructure<dim>::displacement(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorDim::Zero();
    }

    /** 
    * @brief Returns the stress at a given point in the inclusion.
    * @param x The position vector.
    * @param NodeType Type of node.
    * @param ElementType Type of element.
    * @param SimplexDim Type of simplex.
    */

    template<int dim>
    typename InclusionMicrostructure<dim>::MatrixDim InclusionMicrostructure<dim>::stress(const VectorDim& x,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        MatrixDim temp(MatrixDim::Zero());
        for(const auto& inclusion : this->eshelbyInclusions() )
        {
            for(const auto& shift : this->ddBase.periodicShifts)
            {
                temp+=inclusion.second->stress(x+shift);
            }
        }
        return temp;
    }

    /**
    * @brief Returns the mobile concentration of the inclusion.
    * @param VectorDimRef The position vector.
    * @param NodeType Type of node.
    * @param ElementType Type of element.
    * @param SimplexDim Type of simplex.  
    */


    template<int dim>
    typename InclusionMicrostructure<dim>::VectorMSize InclusionMicrostructure<dim>::mobileConcentration(const VectorDim&,const NodeType* const,const ElementType* const,const SimplexDim* const) const
    {
        return VectorMSize::Zero();
    }

    /**
    * @brief Returns the average stress of the inclusion. 
    * @return The average stress of the inclusion (zero matrix).
    */

    template<int dim>
    typename InclusionMicrostructure<dim>::MatrixDim InclusionMicrostructure<dim>::averageStress() const
    {
        return MatrixDim::Zero();
    }

    template struct InclusionMicrostructure<3>;
}
#endif
