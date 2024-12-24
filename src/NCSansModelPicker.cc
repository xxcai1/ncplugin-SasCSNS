#include "NCSansModelPicker.hh"
#include  "NCSansHelper.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/NCString.hh"
#include <vector>
#include <iostream>

bool NCP::SansModelPicker::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::SansIsotropic NCP::SansModelPicker::createFromInfo( const NC::Info& info)
{
  auto iq = SansModelPicker(info);
  return SansIsotropic(iq.getQ(), iq.getI());
}

NCP::SansModelPicker::IqCalType NCP::SansModelPicker::getIqCalType(const NC::Info::CustomSectionData& data) const
{
  for(auto line:data)
  {
    if(line.at(0)=="DirectLoad")
      return IqCalType::kDirectLoad;
    else if(line.at(0)=="HardSphere")
      return IqCalType::kHardSphere;
  }
  return IqCalType::kUndefined;
}

void NCP::SansModelPicker::IqHardSphere(const NC::Info::CustomSectionData& data, double sld, double numden)
{
  //radius
  auto it_r=findCustomLineIter(data, "radius");
  if(it_r->size()!=2)
    NCRYSTAL_THROW2( BadInput,"radius in the @CUSTOM_"<<pluginNameUpperCase()
                   <<" radius field should prove one parameter" );
  double radius(0.);
  if ( !NC::safe_str2dbl( it_r->at(1), radius ) )
    NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                    <<" section radius" );

  double R=radius;
  double R3=R*R*R;
  double V = 4./3.*NC::kPi*R3;
  double atomNumInSphere = V*numden;


  m_Q=NC::logspace(-5, 2,100000);
  m_I.reserve(m_Q.size());

  for(double q:m_Q)
  {
    if(q<1e-5) // approximate by the limit at zero, fixme: this should be found automatically
      m_I.push_back(1.77777777777777777777777778*NC::kPi*NC::kPi*pow(radius,6)*sld*sld/atomNumInSphere);
    else
    {
      double P = 3* (sin(q*R) - q*R*cos(q*R))/(R3* q*q*q);
      m_I.push_back( V*V* sld* sld* P*P/atomNumInSphere);
    }
  }

}

void NCP::SansModelPicker::IqDirectLoad(const NC::Info::CustomSectionData& data)
{
  //Verify we have 3 lines and 2 vectors has identical number of elements
  // if ( data.size() != 2 || data.at(0).size()!=data.at(1).size() )
  //   NCRYSTAL_THROW2(BadInput,"Data in the @CUSTOM_"<<pluginNameUpperCase()
  //                   <<" section should contain two lines that describing an I(Q) (scattering intensity function)");

  auto dataQ=findCustomLineIter(data, "Q");
  auto dataI=findCustomLineIter(data, "I");

  //Parse and validate values:
  m_I.reserve(dataQ->size());
  m_Q.reserve(dataQ->size());

  for(unsigned i=1;i<dataQ->size();i++)
  {
    double temp(0.);
    if ( !NC::safe_str2dbl( dataI->at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section I" );
    m_I.push_back(temp);

    if ( !NC::safe_str2dbl( dataQ->at(i), temp ) )
      NCRYSTAL_THROW2( BadInput,"Invalid values specified in the @CUSTOM_"<<pluginNameUpperCase()
                       <<" section Q" );
    m_Q.push_back(temp);
  }
}

NCP::SansModelPicker::SansModelPicker( const NC::Info& info )
{
  // std::cout << "data source " << info.getDataSourceName() << ", phases " << info.getPhases().size() << std::endl;
  // std::cout << "Info density " << info.getNumberDensity().dbl() << std::endl;

  auto phaselist = info.getPhases();
  const NC::Info *sansinfo = nullptr;
  double sld = 0;
  double numden = 0;
  if(phaselist.empty())
  {
    sld = info.getSLD().dbl() * 0.01 ;  //in sqrt(barn)
    numden = info.getNumberDensity().dbl();  // in atoms/Aa^3
    sansinfo = &info;
  }
  else if(phaselist.size()==2)
  {
    for(const auto& phase : phaselist)
    {
      const auto *ainfo = phase.second.get();
      double frac = phase.first;
      unsigned seccnt = ainfo->countCustomSections(pluginNameUpperCase());

      if(seccnt == 1)
      {
        // the phase for the sans
        sansinfo = ainfo;
        numden = ainfo->getNumberDensity().dbl() * frac;
        sld += ainfo->getSLD().dbl() * 0.01 ; 
        std::cout << "sld seccnt == 1, sld " << sld 
                  << ", frac " << frac 
                  << ", numden " << numden << std::endl;
      }
      else if(seccnt==0)
      {
        // the phase for the solvent
        sld -= ainfo->getSLD().dbl() * 0.01 ; 
        std::cout << "sld seccnt == 0, " << sld 
                  << ", frac " << frac 
                  << ", numden " << ainfo->getNumberDensity().dbl() * frac << std::endl;

      }
      else
      {
        NCRYSTAL_THROW2(BadInput," @CUSTOM_"<<pluginNameUpperCase()<< " contains " <<  seccnt << " section");
      }
    }
  }
  else {
    NCRYSTAL_THROW2(BadInput," info.getPhases()" << " contains " <<  phaselist.size() << " Phases");
  }

  if(!sansinfo)
  {
    NCRYSTAL_THROW2(BadInput," can not find the sans phase");
  }



  printf("sdl %g, number density %g\n",sld, numden);


  // unsigned numSec =  info.countCustomSections( pluginNameUpperCase() ); fixme

  NC::Info::CustomSectionData data = sansinfo->getCustomSection( pluginNameUpperCase(), 0 );

  switch(getIqCalType(data)) {
    case kDirectLoad:
      IqDirectLoad(data);
      break;
    case kHardSphere:
      IqHardSphere(data,std::abs(sld),numden);
      break;
    default :
      NCRYSTAL_THROW2(BadInput," @CUSTOM_"<<pluginNameUpperCase()<< " with undefined load method");
  }
}
