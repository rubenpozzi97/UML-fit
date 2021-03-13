// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME AngDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "interface/PdfRT.h"
#include "interface/PdfWT.h"
#include "interface/DecayRate.h"
#include "interface/PdfSigAng.h"
#include "interface/BoundCheck.h"
#include "interface/Penalty.h"
#include "interface/BoundDist.h"
#include "interface/PdfSigAngMass.h"
#include "interface/PdfSigMass.h"
#include "interface/ShapeSigAng.h"
#include "interface/Fitter.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PdfRT(void *p = 0);
   static void *newArray_PdfRT(Long_t size, void *p);
   static void delete_PdfRT(void *p);
   static void deleteArray_PdfRT(void *p);
   static void destruct_PdfRT(void *p);
   static void streamer_PdfRT(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfRT*)
   {
      ::PdfRT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfRT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfRT", ::PdfRT::Class_Version(), "interface/PdfRT.h", 23,
                  typeid(::PdfRT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfRT::Dictionary, isa_proxy, 16,
                  sizeof(::PdfRT) );
      instance.SetNew(&new_PdfRT);
      instance.SetNewArray(&newArray_PdfRT);
      instance.SetDelete(&delete_PdfRT);
      instance.SetDeleteArray(&deleteArray_PdfRT);
      instance.SetDestructor(&destruct_PdfRT);
      instance.SetStreamerFunc(&streamer_PdfRT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfRT*)
   {
      return GenerateInitInstanceLocal((::PdfRT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfRT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdfWT(void *p = 0);
   static void *newArray_PdfWT(Long_t size, void *p);
   static void delete_PdfWT(void *p);
   static void deleteArray_PdfWT(void *p);
   static void destruct_PdfWT(void *p);
   static void streamer_PdfWT(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfWT*)
   {
      ::PdfWT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfWT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfWT", ::PdfWT::Class_Version(), "interface/PdfWT.h", 23,
                  typeid(::PdfWT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfWT::Dictionary, isa_proxy, 16,
                  sizeof(::PdfWT) );
      instance.SetNew(&new_PdfWT);
      instance.SetNewArray(&newArray_PdfWT);
      instance.SetDelete(&delete_PdfWT);
      instance.SetDeleteArray(&deleteArray_PdfWT);
      instance.SetDestructor(&destruct_PdfWT);
      instance.SetStreamerFunc(&streamer_PdfWT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfWT*)
   {
      return GenerateInitInstanceLocal((::PdfWT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfWT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_DecayRate(void *p = 0);
   static void *newArray_DecayRate(Long_t size, void *p);
   static void delete_DecayRate(void *p);
   static void deleteArray_DecayRate(void *p);
   static void destruct_DecayRate(void *p);
   static void streamer_DecayRate(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DecayRate*)
   {
      ::DecayRate *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DecayRate >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DecayRate", ::DecayRate::Class_Version(), "interface/DecayRate.h", 22,
                  typeid(::DecayRate), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::DecayRate::Dictionary, isa_proxy, 16,
                  sizeof(::DecayRate) );
      instance.SetNew(&new_DecayRate);
      instance.SetNewArray(&newArray_DecayRate);
      instance.SetDelete(&delete_DecayRate);
      instance.SetDeleteArray(&deleteArray_DecayRate);
      instance.SetDestructor(&destruct_DecayRate);
      instance.SetStreamerFunc(&streamer_DecayRate);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DecayRate*)
   {
      return GenerateInitInstanceLocal((::DecayRate*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DecayRate*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdfSigAng(void *p = 0);
   static void *newArray_PdfSigAng(Long_t size, void *p);
   static void delete_PdfSigAng(void *p);
   static void deleteArray_PdfSigAng(void *p);
   static void destruct_PdfSigAng(void *p);
   static void streamer_PdfSigAng(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfSigAng*)
   {
      ::PdfSigAng *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfSigAng >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfSigAng", ::PdfSigAng::Class_Version(), "interface/PdfSigAng.h", 23,
                  typeid(::PdfSigAng), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfSigAng::Dictionary, isa_proxy, 16,
                  sizeof(::PdfSigAng) );
      instance.SetNew(&new_PdfSigAng);
      instance.SetNewArray(&newArray_PdfSigAng);
      instance.SetDelete(&delete_PdfSigAng);
      instance.SetDeleteArray(&deleteArray_PdfSigAng);
      instance.SetDestructor(&destruct_PdfSigAng);
      instance.SetStreamerFunc(&streamer_PdfSigAng);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfSigAng*)
   {
      return GenerateInitInstanceLocal((::PdfSigAng*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfSigAng*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_BoundCheck(void *p = 0);
   static void *newArray_BoundCheck(Long_t size, void *p);
   static void delete_BoundCheck(void *p);
   static void deleteArray_BoundCheck(void *p);
   static void destruct_BoundCheck(void *p);
   static void streamer_BoundCheck(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BoundCheck*)
   {
      ::BoundCheck *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::BoundCheck >(0);
      static ::ROOT::TGenericClassInfo 
         instance("BoundCheck", ::BoundCheck::Class_Version(), "interface/BoundCheck.h", 17,
                  typeid(::BoundCheck), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::BoundCheck::Dictionary, isa_proxy, 16,
                  sizeof(::BoundCheck) );
      instance.SetNew(&new_BoundCheck);
      instance.SetNewArray(&newArray_BoundCheck);
      instance.SetDelete(&delete_BoundCheck);
      instance.SetDeleteArray(&deleteArray_BoundCheck);
      instance.SetDestructor(&destruct_BoundCheck);
      instance.SetStreamerFunc(&streamer_BoundCheck);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BoundCheck*)
   {
      return GenerateInitInstanceLocal((::BoundCheck*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::BoundCheck*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Penalty(void *p = 0);
   static void *newArray_Penalty(Long_t size, void *p);
   static void delete_Penalty(void *p);
   static void deleteArray_Penalty(void *p);
   static void destruct_Penalty(void *p);
   static void streamer_Penalty(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Penalty*)
   {
      ::Penalty *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Penalty >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Penalty", ::Penalty::Class_Version(), "interface/Penalty.h", 17,
                  typeid(::Penalty), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Penalty::Dictionary, isa_proxy, 16,
                  sizeof(::Penalty) );
      instance.SetNew(&new_Penalty);
      instance.SetNewArray(&newArray_Penalty);
      instance.SetDelete(&delete_Penalty);
      instance.SetDeleteArray(&deleteArray_Penalty);
      instance.SetDestructor(&destruct_Penalty);
      instance.SetStreamerFunc(&streamer_Penalty);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Penalty*)
   {
      return GenerateInitInstanceLocal((::Penalty*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Penalty*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_BoundDist(void *p = 0);
   static void *newArray_BoundDist(Long_t size, void *p);
   static void delete_BoundDist(void *p);
   static void deleteArray_BoundDist(void *p);
   static void destruct_BoundDist(void *p);
   static void streamer_BoundDist(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BoundDist*)
   {
      ::BoundDist *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::BoundDist >(0);
      static ::ROOT::TGenericClassInfo 
         instance("BoundDist", ::BoundDist::Class_Version(), "interface/BoundDist.h", 18,
                  typeid(::BoundDist), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::BoundDist::Dictionary, isa_proxy, 16,
                  sizeof(::BoundDist) );
      instance.SetNew(&new_BoundDist);
      instance.SetNewArray(&newArray_BoundDist);
      instance.SetDelete(&delete_BoundDist);
      instance.SetDeleteArray(&deleteArray_BoundDist);
      instance.SetDestructor(&destruct_BoundDist);
      instance.SetStreamerFunc(&streamer_BoundDist);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BoundDist*)
   {
      return GenerateInitInstanceLocal((::BoundDist*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::BoundDist*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdfSigAngMass(void *p = 0);
   static void *newArray_PdfSigAngMass(Long_t size, void *p);
   static void delete_PdfSigAngMass(void *p);
   static void deleteArray_PdfSigAngMass(void *p);
   static void destruct_PdfSigAngMass(void *p);
   static void streamer_PdfSigAngMass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfSigAngMass*)
   {
      ::PdfSigAngMass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfSigAngMass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfSigAngMass", ::PdfSigAngMass::Class_Version(), "interface/PdfSigAngMass.h", 25,
                  typeid(::PdfSigAngMass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfSigAngMass::Dictionary, isa_proxy, 16,
                  sizeof(::PdfSigAngMass) );
      instance.SetNew(&new_PdfSigAngMass);
      instance.SetNewArray(&newArray_PdfSigAngMass);
      instance.SetDelete(&delete_PdfSigAngMass);
      instance.SetDeleteArray(&deleteArray_PdfSigAngMass);
      instance.SetDestructor(&destruct_PdfSigAngMass);
      instance.SetStreamerFunc(&streamer_PdfSigAngMass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfSigAngMass*)
   {
      return GenerateInitInstanceLocal((::PdfSigAngMass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfSigAngMass*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_PdfSigMass(void *p = 0);
   static void *newArray_PdfSigMass(Long_t size, void *p);
   static void delete_PdfSigMass(void *p);
   static void deleteArray_PdfSigMass(void *p);
   static void destruct_PdfSigMass(void *p);
   static void streamer_PdfSigMass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PdfSigMass*)
   {
      ::PdfSigMass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PdfSigMass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PdfSigMass", ::PdfSigMass::Class_Version(), "interface/PdfSigMass.h", 26,
                  typeid(::PdfSigMass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PdfSigMass::Dictionary, isa_proxy, 16,
                  sizeof(::PdfSigMass) );
      instance.SetNew(&new_PdfSigMass);
      instance.SetNewArray(&newArray_PdfSigMass);
      instance.SetDelete(&delete_PdfSigMass);
      instance.SetDeleteArray(&deleteArray_PdfSigMass);
      instance.SetDestructor(&destruct_PdfSigMass);
      instance.SetStreamerFunc(&streamer_PdfSigMass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PdfSigMass*)
   {
      return GenerateInitInstanceLocal((::PdfSigMass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PdfSigMass*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ShapeSigAng(void *p = 0);
   static void *newArray_ShapeSigAng(Long_t size, void *p);
   static void delete_ShapeSigAng(void *p);
   static void deleteArray_ShapeSigAng(void *p);
   static void destruct_ShapeSigAng(void *p);
   static void streamer_ShapeSigAng(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ShapeSigAng*)
   {
      ::ShapeSigAng *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ShapeSigAng >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ShapeSigAng", ::ShapeSigAng::Class_Version(), "interface/ShapeSigAng.h", 23,
                  typeid(::ShapeSigAng), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ShapeSigAng::Dictionary, isa_proxy, 16,
                  sizeof(::ShapeSigAng) );
      instance.SetNew(&new_ShapeSigAng);
      instance.SetNewArray(&newArray_ShapeSigAng);
      instance.SetDelete(&delete_ShapeSigAng);
      instance.SetDeleteArray(&deleteArray_ShapeSigAng);
      instance.SetDestructor(&destruct_ShapeSigAng);
      instance.SetStreamerFunc(&streamer_ShapeSigAng);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ShapeSigAng*)
   {
      return GenerateInitInstanceLocal((::ShapeSigAng*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ShapeSigAng*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Fitter(void *p = 0);
   static void *newArray_Fitter(Long_t size, void *p);
   static void delete_Fitter(void *p);
   static void deleteArray_Fitter(void *p);
   static void destruct_Fitter(void *p);
   static void streamer_Fitter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Fitter*)
   {
      ::Fitter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Fitter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Fitter", ::Fitter::Class_Version(), "interface/Fitter.h", 24,
                  typeid(::Fitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Fitter::Dictionary, isa_proxy, 16,
                  sizeof(::Fitter) );
      instance.SetNew(&new_Fitter);
      instance.SetNewArray(&newArray_Fitter);
      instance.SetDelete(&delete_Fitter);
      instance.SetDeleteArray(&deleteArray_Fitter);
      instance.SetDestructor(&destruct_Fitter);
      instance.SetStreamerFunc(&streamer_Fitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Fitter*)
   {
      return GenerateInitInstanceLocal((::Fitter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Fitter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PdfRT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfRT::Class_Name()
{
   return "PdfRT";
}

//______________________________________________________________________________
const char *PdfRT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfRT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfRT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfRT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfRT*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdfWT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfWT::Class_Name()
{
   return "PdfWT";
}

//______________________________________________________________________________
const char *PdfWT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfWT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfWT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfWT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfWT*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr DecayRate::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *DecayRate::Class_Name()
{
   return "DecayRate";
}

//______________________________________________________________________________
const char *DecayRate::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DecayRate::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *DecayRate::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *DecayRate::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DecayRate*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdfSigAng::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfSigAng::Class_Name()
{
   return "PdfSigAng";
}

//______________________________________________________________________________
const char *PdfSigAng::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfSigAng::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfSigAng::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfSigAng::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAng*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr BoundCheck::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *BoundCheck::Class_Name()
{
   return "BoundCheck";
}

//______________________________________________________________________________
const char *BoundCheck::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BoundCheck*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int BoundCheck::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BoundCheck*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *BoundCheck::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BoundCheck*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *BoundCheck::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BoundCheck*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Penalty::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Penalty::Class_Name()
{
   return "Penalty";
}

//______________________________________________________________________________
const char *Penalty::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Penalty*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Penalty::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Penalty*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Penalty::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Penalty*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Penalty::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Penalty*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr BoundDist::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *BoundDist::Class_Name()
{
   return "BoundDist";
}

//______________________________________________________________________________
const char *BoundDist::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BoundDist*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int BoundDist::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::BoundDist*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *BoundDist::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BoundDist*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *BoundDist::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::BoundDist*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdfSigAngMass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfSigAngMass::Class_Name()
{
   return "PdfSigAngMass";
}

//______________________________________________________________________________
const char *PdfSigAngMass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAngMass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfSigAngMass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAngMass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfSigAngMass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAngMass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfSigAngMass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigAngMass*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr PdfSigMass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PdfSigMass::Class_Name()
{
   return "PdfSigMass";
}

//______________________________________________________________________________
const char *PdfSigMass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigMass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PdfSigMass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PdfSigMass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PdfSigMass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigMass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PdfSigMass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PdfSigMass*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ShapeSigAng::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ShapeSigAng::Class_Name()
{
   return "ShapeSigAng";
}

//______________________________________________________________________________
const char *ShapeSigAng::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ShapeSigAng*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ShapeSigAng::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ShapeSigAng*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ShapeSigAng::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ShapeSigAng*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ShapeSigAng::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ShapeSigAng*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Fitter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Fitter::Class_Name()
{
   return "Fitter";
}

//______________________________________________________________________________
const char *Fitter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Fitter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Fitter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Fitter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Fitter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Fitter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Fitter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Fitter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PdfRT::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfRT.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, PdfRT::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfRT::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfRT(void *p) {
      return  p ? new(p) ::PdfRT : new ::PdfRT;
   }
   static void *newArray_PdfRT(Long_t nElements, void *p) {
      return p ? new(p) ::PdfRT[nElements] : new ::PdfRT[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfRT(void *p) {
      delete ((::PdfRT*)p);
   }
   static void deleteArray_PdfRT(void *p) {
      delete [] ((::PdfRT*)p);
   }
   static void destruct_PdfRT(void *p) {
      typedef ::PdfRT current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfRT(TBuffer &buf, void *obj) {
      ((::PdfRT*)obj)->::PdfRT::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfRT

//______________________________________________________________________________
void PdfWT::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfWT.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, PdfWT::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfWT::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfWT(void *p) {
      return  p ? new(p) ::PdfWT : new ::PdfWT;
   }
   static void *newArray_PdfWT(Long_t nElements, void *p) {
      return p ? new(p) ::PdfWT[nElements] : new ::PdfWT[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfWT(void *p) {
      delete ((::PdfWT*)p);
   }
   static void deleteArray_PdfWT(void *p) {
      delete [] ((::PdfWT*)p);
   }
   static void destruct_PdfWT(void *p) {
      typedef ::PdfWT current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfWT(TBuffer &buf, void *obj) {
      ((::PdfWT*)obj)->::PdfWT::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfWT

//______________________________________________________________________________
void DecayRate::Streamer(TBuffer &R__b)
{
   // Stream an object of class DecayRate.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      PenTerm.Streamer(R__b);
      R__b >> isPenalised;
      R__b.CheckByteCount(R__s, R__c, DecayRate::IsA());
   } else {
      R__c = R__b.WriteVersion(DecayRate::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      PenTerm.Streamer(R__b);
      R__b << isPenalised;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DecayRate(void *p) {
      return  p ? new(p) ::DecayRate : new ::DecayRate;
   }
   static void *newArray_DecayRate(Long_t nElements, void *p) {
      return p ? new(p) ::DecayRate[nElements] : new ::DecayRate[nElements];
   }
   // Wrapper around operator delete
   static void delete_DecayRate(void *p) {
      delete ((::DecayRate*)p);
   }
   static void deleteArray_DecayRate(void *p) {
      delete [] ((::DecayRate*)p);
   }
   static void destruct_DecayRate(void *p) {
      typedef ::DecayRate current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_DecayRate(TBuffer &buf, void *obj) {
      ((::DecayRate*)obj)->::DecayRate::Streamer(buf);
   }
} // end of namespace ROOT for class ::DecayRate

//______________________________________________________________________________
void PdfSigAng::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfSigAng.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      mFrac.Streamer(R__b);
      rtAngTerm.Streamer(R__b);
      wtAngTerm.Streamer(R__b);
      PenTerm.Streamer(R__b);
      R__b >> isPenalised;
      R__b.CheckByteCount(R__s, R__c, PdfSigAng::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfSigAng::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      mFrac.Streamer(R__b);
      rtAngTerm.Streamer(R__b);
      wtAngTerm.Streamer(R__b);
      PenTerm.Streamer(R__b);
      R__b << isPenalised;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfSigAng(void *p) {
      return  p ? new(p) ::PdfSigAng : new ::PdfSigAng;
   }
   static void *newArray_PdfSigAng(Long_t nElements, void *p) {
      return p ? new(p) ::PdfSigAng[nElements] : new ::PdfSigAng[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfSigAng(void *p) {
      delete ((::PdfSigAng*)p);
   }
   static void deleteArray_PdfSigAng(void *p) {
      delete [] ((::PdfSigAng*)p);
   }
   static void destruct_PdfSigAng(void *p) {
      typedef ::PdfSigAng current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfSigAng(TBuffer &buf, void *obj) {
      ((::PdfSigAng*)obj)->::PdfSigAng::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfSigAng

//______________________________________________________________________________
void BoundCheck::Streamer(TBuffer &R__b)
{
   // Stream an object of class BoundCheck.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsReal::Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b >> useCTL4;
      R__b >> useCTL15;
      R__b >> verbose;
      R__b.CheckByteCount(R__s, R__c, BoundCheck::IsA());
   } else {
      R__c = R__b.WriteVersion(BoundCheck::IsA(), kTRUE);
      RooAbsReal::Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b << useCTL4;
      R__b << useCTL15;
      R__b << verbose;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_BoundCheck(void *p) {
      return  p ? new(p) ::BoundCheck : new ::BoundCheck;
   }
   static void *newArray_BoundCheck(Long_t nElements, void *p) {
      return p ? new(p) ::BoundCheck[nElements] : new ::BoundCheck[nElements];
   }
   // Wrapper around operator delete
   static void delete_BoundCheck(void *p) {
      delete ((::BoundCheck*)p);
   }
   static void deleteArray_BoundCheck(void *p) {
      delete [] ((::BoundCheck*)p);
   }
   static void destruct_BoundCheck(void *p) {
      typedef ::BoundCheck current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_BoundCheck(TBuffer &buf, void *obj) {
      ((::BoundCheck*)obj)->::BoundCheck::Streamer(buf);
   }
} // end of namespace ROOT for class ::BoundCheck

//______________________________________________________________________________
void Penalty::Streamer(TBuffer &R__b)
{
   // Stream an object of class Penalty.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsReal::Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b >> coeff1;
      R__b >> coeff4;
      R__b >> coeff5;
      R__b >> power;
      R__b >> verbose;
      R__b.CheckByteCount(R__s, R__c, Penalty::IsA());
   } else {
      R__c = R__b.WriteVersion(Penalty::IsA(), kTRUE);
      RooAbsReal::Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b << coeff1;
      R__b << coeff4;
      R__b << coeff5;
      R__b << power;
      R__b << verbose;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Penalty(void *p) {
      return  p ? new(p) ::Penalty : new ::Penalty;
   }
   static void *newArray_Penalty(Long_t nElements, void *p) {
      return p ? new(p) ::Penalty[nElements] : new ::Penalty[nElements];
   }
   // Wrapper around operator delete
   static void delete_Penalty(void *p) {
      delete ((::Penalty*)p);
   }
   static void deleteArray_Penalty(void *p) {
      delete [] ((::Penalty*)p);
   }
   static void destruct_Penalty(void *p) {
      typedef ::Penalty current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Penalty(TBuffer &buf, void *obj) {
      ((::Penalty*)obj)->::Penalty::Streamer(buf);
   }
} // end of namespace ROOT for class ::Penalty

//______________________________________________________________________________
void BoundDist::Streamer(TBuffer &R__b)
{
   // Stream an object of class BoundDist.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsReal::Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b >> useCTL4;
      R__b >> useCTL15;
      R__b >> verbose;
      R__b.CheckByteCount(R__s, R__c, BoundDist::IsA());
   } else {
      R__c = R__b.WriteVersion(BoundDist::IsA(), kTRUE);
      RooAbsReal::Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      R__b << useCTL4;
      R__b << useCTL15;
      R__b << verbose;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_BoundDist(void *p) {
      return  p ? new(p) ::BoundDist : new ::BoundDist;
   }
   static void *newArray_BoundDist(Long_t nElements, void *p) {
      return p ? new(p) ::BoundDist[nElements] : new ::BoundDist[nElements];
   }
   // Wrapper around operator delete
   static void delete_BoundDist(void *p) {
      delete ((::BoundDist*)p);
   }
   static void deleteArray_BoundDist(void *p) {
      delete [] ((::BoundDist*)p);
   }
   static void destruct_BoundDist(void *p) {
      typedef ::BoundDist current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_BoundDist(TBuffer &buf, void *obj) {
      ((::BoundDist*)obj)->::BoundDist::Streamer(buf);
   }
} // end of namespace ROOT for class ::BoundDist

//______________________________________________________________________________
void PdfSigAngMass::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfSigAngMass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      m.Streamer(R__b);
      mean_rt.Streamer(R__b);
      sigma_rt1.Streamer(R__b);
      sigma_rt2.Streamer(R__b);
      alpha_rt1.Streamer(R__b);
      alpha_rt2.Streamer(R__b);
      n_rt1.Streamer(R__b);
      n_rt2.Streamer(R__b);
      f1rt.Streamer(R__b);
      mean_wt.Streamer(R__b);
      sigma_wt1.Streamer(R__b);
      alpha_wt1.Streamer(R__b);
      alpha_wt2.Streamer(R__b);
      n_wt1.Streamer(R__b);
      n_wt2.Streamer(R__b);
      mFrac.Streamer(R__b);
      constrTerm.Streamer(R__b);
      PenTerm.Streamer(R__b);
      rtAngTerm.Streamer(R__b);
      wtAngTerm.Streamer(R__b);
      rtMassTerm.Streamer(R__b);
      wtMassTerm.Streamer(R__b);
      R__b >> isPenalised;
      R__b.CheckByteCount(R__s, R__c, PdfSigAngMass::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfSigAngMass::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      m.Streamer(R__b);
      mean_rt.Streamer(R__b);
      sigma_rt1.Streamer(R__b);
      sigma_rt2.Streamer(R__b);
      alpha_rt1.Streamer(R__b);
      alpha_rt2.Streamer(R__b);
      n_rt1.Streamer(R__b);
      n_rt2.Streamer(R__b);
      f1rt.Streamer(R__b);
      mean_wt.Streamer(R__b);
      sigma_wt1.Streamer(R__b);
      alpha_wt1.Streamer(R__b);
      alpha_wt2.Streamer(R__b);
      n_wt1.Streamer(R__b);
      n_wt2.Streamer(R__b);
      mFrac.Streamer(R__b);
      constrTerm.Streamer(R__b);
      PenTerm.Streamer(R__b);
      rtAngTerm.Streamer(R__b);
      wtAngTerm.Streamer(R__b);
      rtMassTerm.Streamer(R__b);
      wtMassTerm.Streamer(R__b);
      R__b << isPenalised;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfSigAngMass(void *p) {
      return  p ? new(p) ::PdfSigAngMass : new ::PdfSigAngMass;
   }
   static void *newArray_PdfSigAngMass(Long_t nElements, void *p) {
      return p ? new(p) ::PdfSigAngMass[nElements] : new ::PdfSigAngMass[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfSigAngMass(void *p) {
      delete ((::PdfSigAngMass*)p);
   }
   static void deleteArray_PdfSigAngMass(void *p) {
      delete [] ((::PdfSigAngMass*)p);
   }
   static void destruct_PdfSigAngMass(void *p) {
      typedef ::PdfSigAngMass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfSigAngMass(TBuffer &buf, void *obj) {
      ((::PdfSigAngMass*)obj)->::PdfSigAngMass::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfSigAngMass

//______________________________________________________________________________
void PdfSigMass::Streamer(TBuffer &R__b)
{
   // Stream an object of class PdfSigMass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      mean_rt.Streamer(R__b);
      sigma_rt1.Streamer(R__b);
      sigma_rt2.Streamer(R__b);
      alpha_rt1.Streamer(R__b);
      alpha_rt2.Streamer(R__b);
      n_rt1.Streamer(R__b);
      n_rt2.Streamer(R__b);
      f1rt.Streamer(R__b);
      mean_wt.Streamer(R__b);
      sigma_wt1.Streamer(R__b);
      alpha_wt1.Streamer(R__b);
      alpha_wt2.Streamer(R__b);
      n_wt1.Streamer(R__b);
      n_wt2.Streamer(R__b);
      mFrac.Streamer(R__b);
      rtMassTerm.Streamer(R__b);
      wtMassTerm.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, PdfSigMass::IsA());
   } else {
      R__c = R__b.WriteVersion(PdfSigMass::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      mean_rt.Streamer(R__b);
      sigma_rt1.Streamer(R__b);
      sigma_rt2.Streamer(R__b);
      alpha_rt1.Streamer(R__b);
      alpha_rt2.Streamer(R__b);
      n_rt1.Streamer(R__b);
      n_rt2.Streamer(R__b);
      f1rt.Streamer(R__b);
      mean_wt.Streamer(R__b);
      sigma_wt1.Streamer(R__b);
      alpha_wt1.Streamer(R__b);
      alpha_wt2.Streamer(R__b);
      n_wt1.Streamer(R__b);
      n_wt2.Streamer(R__b);
      mFrac.Streamer(R__b);
      rtMassTerm.Streamer(R__b);
      wtMassTerm.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PdfSigMass(void *p) {
      return  p ? new(p) ::PdfSigMass : new ::PdfSigMass;
   }
   static void *newArray_PdfSigMass(Long_t nElements, void *p) {
      return p ? new(p) ::PdfSigMass[nElements] : new ::PdfSigMass[nElements];
   }
   // Wrapper around operator delete
   static void delete_PdfSigMass(void *p) {
      delete ((::PdfSigMass*)p);
   }
   static void deleteArray_PdfSigMass(void *p) {
      delete [] ((::PdfSigMass*)p);
   }
   static void destruct_PdfSigMass(void *p) {
      typedef ::PdfSigMass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PdfSigMass(TBuffer &buf, void *obj) {
      ((::PdfSigMass*)obj)->::PdfSigMass::Streamer(buf);
   }
} // end of namespace ROOT for class ::PdfSigMass

//______________________________________________________________________________
void ShapeSigAng::Streamer(TBuffer &R__b)
{
   // Stream an object of class ShapeSigAng.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsReal::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b >> isC;
      R__b.CheckByteCount(R__s, R__c, ShapeSigAng::IsA());
   } else {
      R__c = R__b.WriteVersion(ShapeSigAng::IsA(), kTRUE);
      RooAbsReal::Streamer(R__b);
      ctK.Streamer(R__b);
      ctL.Streamer(R__b);
      phi.Streamer(R__b);
      Fl.Streamer(R__b);
      P1.Streamer(R__b);
      P2.Streamer(R__b);
      P3.Streamer(R__b);
      P4p.Streamer(R__b);
      P5p.Streamer(R__b);
      P6p.Streamer(R__b);
      P8p.Streamer(R__b);
      Eff.Streamer(R__b);
      {
         vector<double> &R__stl =  intPart;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b << isC;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ShapeSigAng(void *p) {
      return  p ? new(p) ::ShapeSigAng : new ::ShapeSigAng;
   }
   static void *newArray_ShapeSigAng(Long_t nElements, void *p) {
      return p ? new(p) ::ShapeSigAng[nElements] : new ::ShapeSigAng[nElements];
   }
   // Wrapper around operator delete
   static void delete_ShapeSigAng(void *p) {
      delete ((::ShapeSigAng*)p);
   }
   static void deleteArray_ShapeSigAng(void *p) {
      delete [] ((::ShapeSigAng*)p);
   }
   static void destruct_ShapeSigAng(void *p) {
      typedef ::ShapeSigAng current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ShapeSigAng(TBuffer &buf, void *obj) {
      ((::ShapeSigAng*)obj)->::ShapeSigAng::Streamer(buf);
   }
} // end of namespace ROOT for class ::ShapeSigAng

//______________________________________________________________________________
void Fitter::Streamer(TBuffer &R__b)
{
   // Stream an object of class Fitter.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      name.Streamer(R__b);
      title.Streamer(R__b);
      angPars.Streamer(R__b);
      R__b >> combData;
      R__b >> simPdf;
      R__b >> simPdf_penalty;
      R__b >> boundary;
      R__b >> bound_dist;
      R__b >> penTerm;
      R__b >> constrVars;
      R__b >> fac1;
      R__b >> fac4;
      R__b >> base1;
      R__b >> base4;
      R__b >> max1;
      R__b >> max4;
      R__b >> base1_corr;
      R__b >> base4_corr;
      R__b >> min_base;
      R__b >> maxCoeff;
      R__b >> coeff1;
      R__b >> coeff4;
      R__b >> coeff5;
      R__b >> usedPenalty;
      R__b >> result_free;
      R__b >> result_penalty;
      R__b >> nll;
      R__b >> nll_penalty;
      R__b >> boundDist;
      R__b >> boundDistTime;
      R__b >> minParError;
      R__b >> widthScale;
      {
         vector<Double_t> &R__stl =  vFitResult;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<Double_t> &R__stl =  vFitErrLow;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<Double_t> &R__stl =  vFitErrHigh;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<Double_t> &R__stl =  vImprovResult;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<Double_t> &R__stl =  vResult;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<Double_t> &R__stl =  vConfInterLow;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<Double_t> &R__stl =  vConfInterHigh;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, Fitter::IsA());
   } else {
      R__c = R__b.WriteVersion(Fitter::IsA(), kTRUE);
      name.Streamer(R__b);
      title.Streamer(R__b);
      angPars.Streamer(R__b);
      R__b << combData;
      R__b << simPdf;
      R__b << simPdf_penalty;
      R__b << boundary;
      R__b << bound_dist;
      R__b << penTerm;
      R__b << constrVars;
      R__b << fac1;
      R__b << fac4;
      R__b << base1;
      R__b << base4;
      R__b << max1;
      R__b << max4;
      R__b << base1_corr;
      R__b << base4_corr;
      R__b << min_base;
      R__b << maxCoeff;
      R__b << coeff1;
      R__b << coeff4;
      R__b << coeff5;
      R__b << usedPenalty;
      R__b << result_free;
      R__b << result_penalty;
      R__b << nll;
      R__b << nll_penalty;
      R__b << boundDist;
      R__b << boundDistTime;
      R__b << minParError;
      R__b << widthScale;
      {
         vector<Double_t> &R__stl =  vFitResult;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<Double_t> &R__stl =  vFitErrLow;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<Double_t> &R__stl =  vFitErrHigh;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<Double_t> &R__stl =  vImprovResult;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<Double_t> &R__stl =  vResult;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<Double_t> &R__stl =  vConfInterLow;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<Double_t> &R__stl =  vConfInterHigh;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<Double_t>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Fitter(void *p) {
      return  p ? new(p) ::Fitter : new ::Fitter;
   }
   static void *newArray_Fitter(Long_t nElements, void *p) {
      return p ? new(p) ::Fitter[nElements] : new ::Fitter[nElements];
   }
   // Wrapper around operator delete
   static void delete_Fitter(void *p) {
      delete ((::Fitter*)p);
   }
   static void deleteArray_Fitter(void *p) {
      delete [] ((::Fitter*)p);
   }
   static void destruct_Fitter(void *p) {
      typedef ::Fitter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Fitter(TBuffer &buf, void *obj) {
      ((::Fitter*)obj)->::Fitter::Streamer(buf);
   }
} // end of namespace ROOT for class ::Fitter

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 339,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_AngDict_Impl() {
    static const char* headers[] = {
"interface/PdfRT.h",
"interface/PdfWT.h",
"interface/DecayRate.h",
"interface/PdfSigAng.h",
"interface/BoundCheck.h",
"interface/Penalty.h",
"interface/BoundDist.h",
"interface/PdfSigAngMass.h",
"interface/PdfSigMass.h",
"interface/ShapeSigAng.h",
"interface/Fitter.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc820/lcg/root/6.12.07-ikaegh4/include",
"/afs/cern.ch/user/m/mfaria/public/UML-fit/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "AngDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) description)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfRT.h")))  PdfRT;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) description)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfWT.h")))  PdfWT;
class __attribute__((annotate(R"ATTRDUMP(PDF for angular decay rate description)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/DecayRate.h")))  DecayRate;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfSigAng.h")))  PdfSigAng;
class __attribute__((annotate(R"ATTRDUMP(Step function: 0 inside physical region, 1 outside it)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/BoundCheck.h")))  BoundCheck;
class __attribute__((annotate(R"ATTRDUMP(Penalty term in parameter space)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/Penalty.h")))  Penalty;
class __attribute__((annotate(R"ATTRDUMP(Step function: 0 inside physical region, 1 outside it)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/BoundDist.h")))  BoundDist;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfSigAngMass.h")))  PdfSigAngMass;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/PdfSigMass.h")))  PdfSigMass;
class __attribute__((annotate(R"ATTRDUMP(PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/ShapeSigAng.h")))  ShapeSigAng;
class __attribute__((annotate(R"ATTRDUMP(Code to run the fit and statistical uncertainty)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/Fitter.h")))  Fitter;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "AngDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "interface/PdfRT.h"
#include "interface/PdfWT.h"
#include "interface/DecayRate.h"
#include "interface/PdfSigAng.h"
#include "interface/BoundCheck.h"
#include "interface/Penalty.h"
#include "interface/BoundDist.h"
#include "interface/PdfSigAngMass.h"
#include "interface/PdfSigMass.h"
#include "interface/ShapeSigAng.h"
#include "interface/Fitter.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"BoundCheck", payloadCode, "@",
"BoundDist", payloadCode, "@",
"DecayRate", payloadCode, "@",
"Fitter", payloadCode, "@",
"PdfRT", payloadCode, "@",
"PdfSigAng", payloadCode, "@",
"PdfSigAngMass", payloadCode, "@",
"PdfSigMass", payloadCode, "@",
"PdfWT", payloadCode, "@",
"Penalty", payloadCode, "@",
"ShapeSigAng", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("AngDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_AngDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_AngDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_AngDict() {
  TriggerDictionaryInitialization_AngDict_Impl();
}
