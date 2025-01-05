#include <emscripten/bind.h>
#include <emscripten/val.h>
#include "scs.h"
#include <memory>
#include <string>

using namespace emscripten;

// Helper function to convert JS array to C array
template<typename T>
std::vector<T> convertJSArrayToVector(const val& jsArray) {
    auto length = jsArray["length"].as<unsigned>();
    std::vector<T> result(length);
    val memory = val::module_property("HEAPF64");
    for (unsigned i = 0; i < length; ++i) {
        result[i] = jsArray[i].as<T>();
    }
    return result;
}

template<typename T>
T getValueOrDefault(const val& obj, const std::string& key, T defaultValue) {
    if (obj.hasOwnProperty(key.c_str()) && !obj[key.c_str()].isUndefined() && !obj[key.c_str()].isNull()) {
        return obj[key.c_str()].as<T>();
    }
    return defaultValue;
}

val getArrayOrNull(const val& obj, const std::string& key) {
    if (obj.hasOwnProperty(key.c_str()) && !obj[key.c_str()].isUndefined() && !obj[key.c_str()].isNull()) {
        return obj[key.c_str()];
    }
    return val::null();
}

// Helper class to wrap ScsMatrix for easier JS interaction
class ScsMatrixWrapper {
public:
    ScsMatrixWrapper(val x, val i, val p, int m, int n) {
        matrix = std::make_unique<ScsMatrix>();
        x_data = convertJSArrayToVector<scs_float>(x);
        i_data = convertJSArrayToVector<scs_int>(i);
        p_data = convertJSArrayToVector<scs_int>(p);
        
        matrix->x = x_data.data();
        matrix->i = i_data.data();
        matrix->p = p_data.data();
        matrix->m = m;
        matrix->n = n;
    }

    ScsMatrix* get() { return matrix.get(); }

private:
    std::unique_ptr<ScsMatrix> matrix;
    std::vector<scs_float> x_data;
    std::vector<scs_int> i_data;
    std::vector<scs_int> p_data;
};

// Helper class to wrap ScsData
class ScsDataWrapper {
public:
    ScsDataWrapper(int m, int n, val A_x, val A_i, val A_p, 
                  const val& dataObj) {
        data = std::make_unique<ScsData>();
        data->m = m;
        data->n = n;

        // Set up A matrix
        A = std::make_unique<ScsMatrixWrapper>(A_x, A_i, A_p, m, n);
        data->A = A->get();

        // Set up P matrix if provided
        val P_x = getArrayOrNull(dataObj, "P_x");
        val P_i = getArrayOrNull(dataObj, "P_i");
        val P_p = getArrayOrNull(dataObj, "P_p");
        
        if (!P_x.isNull() && !P_i.isNull() && !P_p.isNull()) {
            P = std::make_unique<ScsMatrixWrapper>(P_x, P_i, P_p, n, n);
            data->P = P->get();
        } else {
            data->P = nullptr;
        }

        // Convert b and c vectors
        b_data = convertJSArrayToVector<scs_float>(dataObj["b"]);
        c_data = convertJSArrayToVector<scs_float>(dataObj["c"]);
        data->b = b_data.data();
        data->c = c_data.data();
    }

    ScsData* get() { return data.get(); }

private:
    std::unique_ptr<ScsData> data;
    std::unique_ptr<ScsMatrixWrapper> A;
    std::unique_ptr<ScsMatrixWrapper> P;
    std::vector<scs_float> b_data;
    std::vector<scs_float> c_data;
};

// Helper class to wrap ScsCone
class ScsConeWrapper {
public:
    ScsConeWrapper(const val& coneObj) {
        cone = std::make_unique<ScsCone>();
        
        // Set default values
        cone->z = getValueOrDefault(coneObj, "z", 0);
        cone->l = getValueOrDefault(coneObj, "l", 0);
        cone->bsize = getValueOrDefault(coneObj, "bsize", 0);
        cone->qsize = getValueOrDefault(coneObj, "qsize", 0);
        cone->ssize = getValueOrDefault(coneObj, "ssize", 0);
        cone->ep = getValueOrDefault(coneObj, "ep", 0);
        cone->ed = getValueOrDefault(coneObj, "ed", 0);
        cone->psize = getValueOrDefault(coneObj, "psize", 0);

        // Handle arrays
        val bu = getArrayOrNull(coneObj, "bu");
        val bl = getArrayOrNull(coneObj, "bl");
        if (!bu.isNull() && !bl.isNull()) {
            bu_data = convertJSArrayToVector<scs_float>(bu);
            bl_data = convertJSArrayToVector<scs_float>(bl);
            cone->bu = bu_data.data();
            cone->bl = bl_data.data();
        } else {
            cone->bu = nullptr;
            cone->bl = nullptr;
        }

        val q = getArrayOrNull(coneObj, "q");
        if (!q.isNull()) {
            q_data = convertJSArrayToVector<scs_int>(q);
            cone->q = q_data.data();
        } else {
            cone->q = nullptr;
        }

        val s = getArrayOrNull(coneObj, "s");
        if (!s.isNull()) {
            s_data = convertJSArrayToVector<scs_int>(s);
            cone->s = s_data.data();
        } else {
            cone->s = nullptr;
        }

        val p = getArrayOrNull(coneObj, "p");
        if (!p.isNull()) {
            p_data = convertJSArrayToVector<scs_float>(p);
            cone->p = p_data.data();
        } else {
            cone->p = nullptr;
        }
    }

    ScsCone* get() { return cone.get(); }

private:
    std::unique_ptr<ScsCone> cone;
    std::vector<scs_float> bu_data, bl_data;
    std::vector<scs_int> q_data, s_data;
    std::vector<scs_float> p_data;
};

// Wrapper function for scs_solve that returns solution as a JavaScript object
val solve_scs(ScsDataWrapper& data, ScsConeWrapper& cone, const ScsSettings& settings, val solutionObj = val::null()) {
    auto sol = std::make_unique<ScsSolution>();
    auto info = std::make_unique<ScsInfo>();
    
    std::vector<scs_float> x_data, y_data, s_data;
    if (!solutionObj.isNull()) {
        // Initialize sol with provided solution object for warm starting
        x_data = convertJSArrayToVector<scs_float>(solutionObj["x"]);
        y_data = convertJSArrayToVector<scs_float>(solutionObj["y"]);
        s_data = convertJSArrayToVector<scs_float>(solutionObj["s"]);
        sol->x = x_data.data();
        sol->y = y_data.data();
        sol->s = s_data.data();
    } else {
        // Initialize solution vectors to nullptr
        sol->x = nullptr;
        sol->y = nullptr;
        sol->s = nullptr;
    }

    // Call SCS solver
    scs_int result = scs(data.get(), cone.get(), &settings, sol.get(), info.get());

    // Convert results to JavaScript object
    val solution = val::object();
    if (result >= 0) {  // If solved successfully
        int n = data.get()->n;
        int m = data.get()->m;
        
        // Convert solution vectors to JavaScript arrays
        val x_array = val::array();
        val y_array = val::array();
        val s_array = val::array();
        
        for (int i = 0; i < n; i++) x_array.set(i, sol->x[i]);
        for (int i = 0; i < m; i++) y_array.set(i, sol->y[i]);
        for (int i = 0; i < m; i++) s_array.set(i, sol->s[i]);
        
        solution.set("x", x_array);
        solution.set("y", y_array);
        solution.set("s", s_array);
    }
    
    // Add info to result
    val info_obj = val::object();
    info_obj.set("iter", info->iter);
    info_obj.set("pobj", info->pobj);
    info_obj.set("dobj", info->dobj);
    info_obj.set("resPri", info->res_pri);
    info_obj.set("resDual", info->res_dual);
    info_obj.set("resInfeas", info->res_infeas);
    info_obj.set("resUnbdd", info->res_unbdd_a);
    info_obj.set("solveTime", info->solve_time);
    info_obj.set("setupTime", info->setup_time);
    
    solution.set("info", info_obj);
    solution.set("status", result);
    
    return solution;
}

// Wrapper function for scs_solve that accepts JS objects for data and cone
val solve_scs(val dataObj, val coneObj, const ScsSettings& settings, val solutionObj = val::null()) {
    // Convert JS data object to ScsDataWrapper
    ScsDataWrapper data(dataObj["m"].as<int>(), dataObj["n"].as<int>(),
                       dataObj["A_x"], dataObj["A_i"], dataObj["A_p"],
                       dataObj);

    // Use the cone constructor
    ScsConeWrapper cone(coneObj);

    return solve_scs(data, cone, settings, solutionObj);
}

// Overload for scs_solve that does not require a solution object
val solve_scs(val dataObj, val coneObj, const ScsSettings& settings) {
    return solve_scs(dataObj, coneObj, settings, val::null());
}

// Wrapper function for scs_set_default_settings
void set_default_settings(ScsSettings& settings) {
    scs_set_default_settings(&settings);
}

EMSCRIPTEN_BINDINGS(scs_module) {
    class_<ScsSettings>("ScsSettings")
        .constructor<>()
        .property("normalize", &ScsSettings::normalize)
        .property("scale", &ScsSettings::scale)
        .property("adaptiveScale", &ScsSettings::adaptive_scale)
        .property("rhoX", &ScsSettings::rho_x)
        .property("maxIters", &ScsSettings::max_iters)
        .property("epsAbs", &ScsSettings::eps_abs)
        .property("epsRel", &ScsSettings::eps_rel)
        .property("epsInfeas", &ScsSettings::eps_infeas)
        .property("alpha", &ScsSettings::alpha)
        .property("timeLimitSecs", &ScsSettings::time_limit_secs)
        .property("verbose", &ScsSettings::verbose)
        .property("warmStart", &ScsSettings::warm_start)
        ;
    
    class_<ScsMatrixWrapper>("ScsMatrix")
        .constructor<val, val, val, int, int>()
        ;
    
    class_<ScsDataWrapper>("ScsData")
        .constructor<int, int, val, val, val, val>()
        ;
    
    class_<ScsConeWrapper>("ScsCone")
        .constructor<val>()
        ;
    
    function("solve", select_overload<val(val, val, const ScsSettings&, val)>(&solve_scs));
    function("solve", select_overload<val(val, val, const ScsSettings&)>(&solve_scs));
    function("setDefaultSettings", &set_default_settings);
    
}
