/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * Written in 2024 by Matthew Maron <mlm335@miami.edu>
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityPy_cpp_
#define _model_DislocationMobilityPy_cpp_

#include <DislocationMobilityPy.h>
#include <filesystem>

namespace model
{
#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    DislocationMobilityPy::DislocationMobilityPy(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in) :
    /* init */ DislocationMobilityBase("Using Py mobility for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,mu_SI(material.mu_SI)
    /* init */,Tm(material.Tm)
    /* init */,cs(material.cs_SI)
    /* init */,pyModuleName(pyModuleName_in)
    {// Set up pyModule
        std::filesystem::path modulePath(pyModuleName);
        const std::string moduleStem(modulePath.stem().string());
        const std::string moduleDir(modulePath.parent_path().string());
        std::cout<<"moduleStem="<<moduleStem<<std::endl;
        std::cout<<"moduleDir="<<moduleDir<<std::endl;
        pybind11::module sys = pybind11::module::import("sys");
        pybind11::list path = sys.attr("path");
        path.append(moduleDir);
        pyModule = pybind11::module::import(moduleStem.c_str());
    }

    double DislocationMobilityPy::velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& )
    {
        Eigen::MatrixXd stress(S*mu_SI*1e-9); // GPa
        Eigen::VectorXd burgers(b); // Vector Direction
        Eigen::VectorXd tangent(xi); // Vector Direction
        Eigen::VectorXd normal(n); // Vector Direction
        const double temp(T); // K
        
        try
        {
           pybind11::object mobilitySolver = pyModule.attr("MobilitySolver")();
           double qpPythonVelocity = mobilitySolver.attr("velocityPy")(stress, burgers, tangent, normal, temp, dL, dt).cast<double>();
           return qpPythonVelocity*100/cs; //Convert from A/ps to code units
        }
        catch (const pybind11::error_already_set& e)
        {
           std::cerr << "Python error: " << e.what() << std::endl;
           std::terminate();
        }
        

    }

#else // COMPILED WITHOUT PYBIND11
    DislocationMobilityPy::DislocationMobilityPy(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in) :
    /* init */ DislocationMobilityBase("Py mobility for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,mu_SI(material.mu_SI)
    /* init */,Tm(material.Tm)
    /* init */,cs(material.cs_SI)
    /* init */,pyModuleName(pyModuleName_in)
    {
    }

    double DislocationMobilityPy::velocity(const MatrixDim& ,
                    const VectorDim& ,
                    const VectorDim& ,
                    const VectorDim& ,
                    const double& ,
                    const double& ,
                    const double& ,
                    const std::shared_ptr<StochasticForceGenerator>& )
    {
        throw std::runtime_error("DislocationMobilityPy used without pybind11");
        return 0.0;
    }
#endif

}
#endif
