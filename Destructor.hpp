#pragma once

// This is a simple class to support calling a function when the object goes
// out of scope. This allows easy RAII cleanup of objects that aren't directly
// supported by std::shared_ptr and std::unique_ptr such as FILE* and file
// descriptors.
struct Destructor
{
    typedef void (*FuncType)();
    FuncType F;
    inline Destructor(FuncType F) : F(F) {}
    inline ~Destructor() { F(); }
};
