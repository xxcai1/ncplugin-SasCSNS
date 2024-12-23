#include "NCSansHelper.hh"


NC::Info::CustomSectionData::const_iterator NCPluginNamespace::findCustomLineIter(const NC::Info::CustomSectionData& data, const std::string& keyword, bool check)
{
  NC::Info::CustomSectionData::const_iterator it(data.end());
  for(auto line=data.begin();line!=data.end();++line)
  {
    if(line->at(0)==keyword)
    {
      it=line;
      break;
    }
  }
  if(check && it==data.end())
    NCRYSTAL_THROW2(BadInput,"findCustomLineIter can not find parameter " << keyword);
  return it;
}
