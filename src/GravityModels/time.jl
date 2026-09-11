## Description #############################################################################
#
# Functions to convert the time representations used in the gravity models.
#
############################################################################################

const _DT_J2000 = DateTime(2000, 1, 1, 12, 0, 0)

"""
    _to_j2000_seconds(time::Number) -> Number
    _to_j2000_seconds(time::DateTime) -> Float64

Convert `time` to the number of elapsed seconds [s] from the J2000.0 epoch
(2000-01-01T12:00:00). If `time` is a number, it is assumed to be already expressed in this
representation and it is returned unchanged.
"""
_to_j2000_seconds(time::Number) = time
_to_j2000_seconds(time::DateTime) = Dates.value(time - _DT_J2000) / 1000

"""
    _from_j2000_seconds(time::Number) -> DateTime

Convert `time`, expressed as the number of elapsed seconds [s] from the J2000.0 epoch
(2000-01-01T12:00:00), to a `DateTime` object with millisecond resolution.
"""
function _from_j2000_seconds(time::Number)
    return _DT_J2000 + Dates.Millisecond(round(Int64, 1000 * time))
end
