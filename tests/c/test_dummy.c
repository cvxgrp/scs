#include "test_dummy.h"

bool test_dummy_method(char** str) {
    if (assertEqualsInt(1, 1)) {
        SUCCEED(str);
    } else {
        FAIL_WITH_MESSAGE(str, "failure");
    }
}