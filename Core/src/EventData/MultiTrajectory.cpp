#include "Acts/EventData/MultiTrajectory.hpp"

template class Acts::detail_lt::TrackStateProxy<
    Acts::MultiTrajectory::MeasurementSizeMax, true>;

template class Acts::detail_lt::TrackStateProxy<
    Acts::MultiTrajectory::MeasurementSizeMax, false>;