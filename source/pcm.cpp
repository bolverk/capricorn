#include "pcm.hpp"
#include "utilities.hpp"

PCM::PCM(void) {}

vector<pair<Primitive,Primitive> > PCM::operator()
(const HydroSnapshot& hs) const
{
  class Interpolator: public Index2Member<pair<Primitive, Primitive> >
  {
  public:

    Interpolator(const HydroSnapshot& hs_i):
    hs_(hs_i) {}

    size_t getLength(void) const
    {
      return hs_.grid.size() - 2;
    }

    pair<Primitive,Primitive> operator()(size_t i) const
    {
      return pair<Primitive,Primitive>
	(hs_.cells[i],hs_.cells[i+1]);
    }

  private:
    const HydroSnapshot& hs_;
  } interpolator(hs);

  return serial_generate(interpolator);
}
