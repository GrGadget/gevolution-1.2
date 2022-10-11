#pragma once


namespace gevolution
{
/*
    From Stroustrup's The C++ Programming Language 4ed
    section 13.3.1
*/
template <class F>
struct Final_action
// A final action can be used to add C++ RAII safety to C API.
// For example: I am using a resource of type Res, that needs to be initialized and released after
// use with a C api, that would be
// my_function()
// {
//     Res X;
//     X = Res_initialize(...);
//     // ...
//     Res_free(X);    
// }
// To achive automatic free of the resource to yield exception safe code one could write
// my_function()
// {
//     Res X;
//     X = Res_initialize(...);
//     auto cleanup = finally([&](){
//         Res_free(X);    
//     });
//     // ...
// }
// At destruction, the object cleanup will execute the lambda passed to the factory.
{
    F clean;
    Final_action(F f) : clean{f} {}
    ~Final_action() { clean(); }
};


template <class F>
Final_action<F> finally(F f)
// Final action factory
{
    return Final_action<F>(f);
}

}
