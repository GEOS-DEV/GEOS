//FUNCTION_BEGIN_PARSE
virtual_int
CohesiveZoneBase::UpdateCohesiveZone( const localIndex index,
                                      const R1Tensor& gap,
                                      const R1Tensor& N,
                                      const std::pair< ElementRegionT*, localIndex >& elem0,
                                      const std::pair< ElementRegionT*, localIndex >& elem1,
                                      R1Tensor& traction,
                                      R2Tensor& stiffness )
{
  return 1;
}
