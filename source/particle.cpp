/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2013 Tobias Pfaff, Nils Thuerey
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Particle data functionality
 *
 ******************************************************************************/

#include <fstream>
#include <cstring>
#if NO_ZLIB!=1
#include <zlib.h>
#endif
#include "particle.h"
#include "levelset.h"
#include "mantaio.h"
#include "vortexpart.h"
#include "turbulencepart.h"

using namespace std;
namespace Manta {


ParticleBase::ParticleBase(FluidSolver* parent)
	: PbClass(parent), mAllowCompress(true), mFreePdata(false)
{
}

ParticleBase::~ParticleBase()
{
	// make sure data fields now parent system is deleted
	for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i)
		mPartData[i]->setParticleSys(NULL);

	if(mFreePdata) {
		for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i)
			delete mPartData[i];
	}

}

std::string ParticleBase::infoString() const
{
	return "ParticleSystem " + mName + " <no info>";
}

void ParticleBase::cloneParticleData(ParticleBase* nm)
{
	// clone additional data , and make sure the copied particle system deletes it
	nm->mFreePdata = true;
	for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i) {
		ParticleDataBase* pdata = mPartData[i]->clone();
		nm->registerPdata(pdata);
	}
}

void ParticleBase::deregister(ParticleDataBase* pdata)
{
	bool done = false;
	// remove pointer from particle data list
	for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i) {
		if(mPartData[i] == pdata) {
			if(i<(IndexInt)mPartData.size()-1)
				mPartData[i] = mPartData[mPartData.size()-1];
			mPartData.pop_back();
			done = true;
		}
	}
	if(!done)
		errMsg("Invalid pointer given, not registered!");
}

// create and attach a new pdata field to this particle system
PbClass* ParticleBase::create(PbType t, PbTypeVec T, const string& name)
{
#	if NOPYTHON!=1
	_args.add("nocheck",true);
	if(t.str() == "") errMsg("Specify particle data type to create");
	//debMsg( "Pdata creating '"<< t.str <<" with size "<< this->getSizeSlow(), 5 );

	PbClass* pyObj = PbClass::createPyObject(t.str() + T.str(), name, _args, this->getParent());

	ParticleDataBase* pdata = dynamic_cast<ParticleDataBase*>(pyObj);
	if(!pdata) {
		errMsg("Unable to get particle data pointer from newly created object. Only create ParticleData type with a ParticleSys.creat() call, eg, PdataReal, PdataVec3 etc.");
		delete pyObj;
		return NULL;
	} else {
		this->registerPdata(pdata);
	}

	// directly init size of new pdata field:
	pdata->resize(this->getSizeSlow());
#	else
	PbClass* pyObj = NULL;
#	endif
	return pyObj;
}

void ParticleBase::registerPdata(ParticleDataBase* pdata)
{
	pdata->setParticleSys(this);
	mPartData.push_back(pdata);

	if(pdata->getType() == ParticleDataBase::TypeReal) {
		ParticleDataImpl<Real>* pd = dynamic_cast< ParticleDataImpl<Real>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as real!");
		this->registerPdataReal(pd);
	} else if(pdata->getType() == ParticleDataBase::TypeInt) {
		ParticleDataImpl<int>* pd = dynamic_cast< ParticleDataImpl<int>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as int!");
		this->registerPdataInt(pd);
	} else if(pdata->getType() == ParticleDataBase::TypeVec3) {
		ParticleDataImpl<Vec3>* pd = dynamic_cast< ParticleDataImpl<Vec3>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as vec3!");
		this->registerPdataVec3(pd);
	}
}
void ParticleBase::registerPdataReal(ParticleDataImpl<Real>* pd) { mPdataReal.push_back(pd); }
void ParticleBase::registerPdataVec3(ParticleDataImpl<Vec3>* pd) { mPdataVec3.push_back(pd); }
void ParticleBase::registerPdataInt (ParticleDataImpl<int >* pd) { mPdataInt .push_back(pd); }

void ParticleBase::addAllPdata()
{
	for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i) {
		mPartData[i]->addEntry();
	}
}


BasicParticleSystem::BasicParticleSystem(FluidSolver* parent)
	: ParticleSystem<BasicParticleData>(parent)
{
	this->mAllowCompress = false;
}

// file io

void BasicParticleSystem::writeParticlesText(const string name) const
{
	ofstream ofs(name.c_str());
	if(!ofs.good()) errMsg("can't open file!");
	ofs << this->size()<<", pdata: "<< mPartData.size()<<" ("<<mPdataInt.size()<<","<<mPdataReal.size()<<","<<mPdataVec3.size()<<") \n";
	for(IndexInt i=0; i<this->size(); ++i) {
		ofs << i<<": "<< this->getPos(i) <<" , "<< this->getStatus(i) <<". ";
		for(IndexInt pd=0; pd<(IndexInt)mPdataInt.size() ; ++pd) ofs << mPdataInt [pd]->get(i)<<" ";
		for(IndexInt pd=0; pd<(IndexInt)mPdataReal.size(); ++pd) ofs << mPdataReal[pd]->get(i)<<" ";
		for(IndexInt pd=0; pd<(IndexInt)mPdataVec3.size(); ++pd) ofs << mPdataVec3[pd]->get(i)<<" ";
		ofs << "\n";
	}
	ofs.close();
}

void BasicParticleSystem::writeParticlesRawPositionsGz(const string name) const
{
#	if NO_ZLIB!=1
	gzFile gzf = (gzFile) safeGzopen(name.c_str(), "wb1");
	if(!gzf) errMsg("can't open file "<<name);
	for(IndexInt i=0; i<this->size(); ++i) {
		Vector3D<float> p = toVec3f(this->getPos(i));
		gzwrite(gzf, &p, sizeof(float)*3);
	}
	gzclose(gzf);
#	else
	cout << "file format not supported without zlib" << endl;
#	endif
}

void BasicParticleSystem::writeParticlesRawVelocityGz(const string name) const
{
#	if NO_ZLIB!=1
	gzFile gzf = (gzFile) safeGzopen(name.c_str(), "wb1");
	if (!gzf) errMsg("can't open file "<<name);
	if( mPdataVec3.size() < 1 ) errMsg("no vec3 particle data channel found!");
	// note , assuming particle data vec3 0 is velocity! make optional...
	for(IndexInt i=0; i<this->size(); ++i) {
		Vector3D<float> p = toVec3f(mPdataVec3[0]->get(i));
		gzwrite(gzf, &p, sizeof(float)*3);
	}
	gzclose(gzf);
#	else
	cout << "file format not supported without zlib" << endl;
#	endif
}


int BasicParticleSystem::load(const string name)
{
	if(name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if(ext == ".uni")
		return readParticlesUni(name, this );
	else if (ext == ".vdb") {
		std::vector<PbClass*> parts;
		parts.push_back(this);
		return readObjectsVDB(name, &parts);
	} else if(ext == ".raw") // raw = uni for now
		return readParticlesUni(name, this );
	else
		errMsg("particle '" + name +"' filetype not supported for loading");
	return 0;
}

int BasicParticleSystem::save(const string name)
{
	if(name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if(ext == ".txt")
		this->writeParticlesText(name);
	else if(ext == ".uni")
		return writeParticlesUni(name, this);
	else if(ext == ".raw") // raw = uni for now
		return writeParticlesUni(name, this);
	else if (ext == ".vdb") {
		std::vector<PbClass*> parts;
		parts.push_back(this);
		return writeObjectsVDB(name, &parts);
	// raw data formats, very basic for simple data transfer to other programs
	} else if(ext == ".posgz")
		this->writeParticlesRawPositionsGz(name);
	else if(ext == ".velgz")
		this->writeParticlesRawVelocityGz(name);
	else
		errMsg("particle '" + name +"' filetype not supported for saving");
	return 0;
}

void BasicParticleSystem::printParts(IndexInt start, IndexInt stop, bool printIndex)
{
	std::ostringstream sstr;
	IndexInt s = (start>0 ? start : 0                      );
	IndexInt e = (stop>0  ? stop  : (IndexInt)mData.size() );
	s = Manta::clamp(s, (IndexInt)0, (IndexInt)mData.size());
	e = Manta::clamp(e, (IndexInt)0, (IndexInt)mData.size());

	for(IndexInt i=s; i<e; ++i) {
		if(printIndex) sstr << i<<": ";
		sstr<<mData[i].pos<<" "<<mData[i].flag<<"\n";
	}
	debMsg( sstr.str() , 1 );
}

std::string BasicParticleSystem::getDataPointer() {
	std::ostringstream out;
	out << &mData;
	return out.str();
}

void BasicParticleSystem::readParticles(BasicParticleSystem* from) {
	// re-allocate all data
	this->resizeAll(from->size());
	assertMsg(from->size() == this->size() , "particle size doesn't match");

	for(int i=0; i<this->size(); ++i) {
		(*this)[i].pos  = (*from)[i].pos;
		(*this)[i].flag = (*from)[i].flag;
	}
	this->transformPositions(from->getParent()->getGridSize(), this->getParent()->getGridSize());
}


// particle data

ParticleDataBase::ParticleDataBase(FluidSolver* parent)
	: PbClass(parent), mpParticleSys(NULL)
{
}

ParticleDataBase::~ParticleDataBase()
{
	// notify parent of deletion
	if(mpParticleSys)
		mpParticleSys->deregister(this);
}


// actual data implementation

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent)
	: ParticleDataBase(parent), mpGridSource(NULL), mGridSourceMAC(false)
{
}

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other)
	: ParticleDataBase(parent), mpGridSource(NULL), mGridSourceMAC(false)
{
	this->mData = other->mData;
	setName(other->getName());
}

template<class T>
ParticleDataImpl<T>::~ParticleDataImpl()
{
}

template<class T>
IndexInt ParticleDataImpl<T>::getSizeSlow() const
{
	return mData.size();
}
template<class T>
void ParticleDataImpl<T>::addEntry()
{
	// add zero'ed entry
	T tmp = T(0.);
	// for debugging, force init:
	//tmp = T(0.02 * mData.size()); // increasing
	//tmp = T(1.); // constant 1
	return mData.push_back(tmp);
}
template<class T>
void ParticleDataImpl<T>::resize(IndexInt s)
{
	mData.resize(s);
}
template<class T>
void ParticleDataImpl<T>::copyValueSlow(IndexInt from, IndexInt to)
{
	this->copyValue(from,to);
}
template<class T>
ParticleDataBase* ParticleDataImpl<T>::clone()
{
	ParticleDataImpl<T>* npd = new ParticleDataImpl<T>(getParent(), this);
	return npd;
}

template<class T>
void ParticleDataImpl<T>::setSource(Grid<T>* grid, bool isMAC)
{
	mpGridSource = grid;
	mGridSourceMAC = isMAC;
	if(isMAC) assertMsg(dynamic_cast<MACGrid*>(grid) != NULL , "Given grid is not a valid MAC grid");
}

template<class T>
void ParticleDataImpl<T>::initNewValue(IndexInt idx, Vec3 pos)
{
	if(!mpGridSource)
		mData[idx] = 0;
	else {
		mData[idx] = mpGridSource->getInterpolated(pos);
	}
}
// special handling needed for velocities
template<>
void ParticleDataImpl<Vec3>::initNewValue(IndexInt idx, Vec3 pos)
{
	if(!mpGridSource)
		mData[idx] = 0;
	else {
		if(!mGridSourceMAC)
			mData[idx] = mpGridSource->getInterpolated(pos);
		else
			mData[idx] = ((MACGrid*)mpGridSource)->getInterpolated(pos);
	}
}

template<typename T>
int ParticleDataImpl<T>::load(string name)
{
	if(name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if(ext == ".uni")
		return readPdataUni<T>(name, this);
	else if (ext == ".vdb") {
		std::vector<PbClass*> parts;
		parts.push_back(this);
		return readObjectsVDB(name, &parts);
	}
	else if(ext == ".raw") // raw = uni for now
		return readPdataUni<T>(name, this);
	else
		errMsg("particle data '" + name +"' filetype not supported for loading");
	return 0;
}

template<typename T>
int ParticleDataImpl<T>::save(string name)
{
	if(name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if(ext == ".uni")
		return writePdataUni<T>(name, this);
	else if (ext == ".vdb") {
		std::vector<PbClass*> parts;
		parts.push_back(this);
		return writeObjectsVDB(name, &parts);
	}
	else if(ext == ".raw") // raw = uni for now
		return writePdataUni<T>(name, this);
	else
		errMsg("particle data '" + name +"' filetype not supported for saving");
	return 0;
}

// specializations

template<>
ParticleDataBase::PdataType ParticleDataImpl<Real>::getType() const
{
	return ParticleDataBase::TypeReal;
}
template<>
ParticleDataBase::PdataType ParticleDataImpl<int>::getType() const
{
	return ParticleDataBase::TypeInt;
}
template<>
ParticleDataBase::PdataType ParticleDataImpl<Vec3>::getType() const
{
	return ParticleDataBase::TypeVec3;
}

// note, we need a flag value for functions such as advection
// ideally, this value should never be modified
int ParticleIndexData::flag = 0;
Vec3 ParticleIndexData::pos = Vec3(0.,0.,0.);

KERNEL(pts) template<class T, class S> void knPdataAdd    (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] += other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataSub    (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] -= other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataMult   (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] *= other[idx]; }
KERNEL(pts) template<class T, class S> void knPdataDiv    (ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other) { me[idx] /= other[idx]; }
KERNEL(pts) template<class T>          void knPdataSafeDiv(ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other) { me[idx] = safeDivide(me[idx], other[idx]); }

KERNEL(pts) template<class T, class S> void knPdataSetScalar    (ParticleDataImpl<T>& me, const S& other) { me[idx]  = other; }
KERNEL(pts) template<class T, class S> void knPdataAddScalar    (ParticleDataImpl<T>& me, const S& other) { me[idx] += other; }
KERNEL(pts) template<class T, class S> void knPdataMultScalar   (ParticleDataImpl<T>& me, const S& other) { me[idx] *= other; }
KERNEL(pts) template<class T, class S> void knPdataScaledAdd    (ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other, const S& factor) { me[idx] += factor * other[idx]; }

KERNEL(pts) template<class T> void knPdataClamp(ParticleDataImpl<T>& me, const T vmin, const T vmax) { me[idx] = clamp(me[idx], vmin, vmax); }
KERNEL(pts) template<class T> void knPdataClampMin(ParticleDataImpl<T>& me, const T vmin)            { me[idx] = std::max(vmin, me[idx]); }
KERNEL(pts) template<class T> void knPdataClampMax(ParticleDataImpl<T>& me, const T vmax)            { me[idx] = std::min(vmax, me[idx]); }
KERNEL(pts)                   void knPdataClampMinVec3(ParticleDataImpl<Vec3>& me, const Real vmin)
{
	me[idx].x = std::max(vmin, me[idx].x);
	me[idx].y = std::max(vmin, me[idx].y);
	me[idx].z = std::max(vmin, me[idx].z);
}
KERNEL(pts)                   void knPdataClampMaxVec3(ParticleDataImpl<Vec3>& me, const Real vmax)
{
	me[idx].x = std::min(vmax, me[idx].x);
	me[idx].y = std::min(vmax, me[idx].y);
	me[idx].z = std::min(vmax, me[idx].z);
}

// python operators


template<typename T>
ParticleDataImpl<T>& ParticleDataImpl<T>::copyFrom(const ParticleDataImpl<T>& a)
{
	assertMsg(a.mData.size() == mData.size() , "different pdata size "<<a.mData.size()<<" vs "<<this->mData.size());
	mData = a.mData;
	return *this;
}

template<typename T>
void ParticleDataImpl<T>::setConst(const T &s)
{
	knPdataSetScalar<T,T> op( *this, s );
}

template<typename T>
void ParticleDataImpl<T>::setConstRange(const T &s, const int begin, const int end)
{
	for(int i=begin; i<end; ++i) (*this)[i] = s;
}

// special set by flag
KERNEL(pts) template<class T, class S> void knPdataSetScalarIntFlag(ParticleDataImpl<T> &me, const S &other, const ParticleDataImpl<int> &t, const int itype)
{
	if(t[idx]&itype) me[idx] = other;
}
template<typename T>
void ParticleDataImpl<T>::setConstIntFlag(const T &s, const ParticleDataImpl<int> &t, const int itype)
{
	knPdataSetScalarIntFlag<T,T> op(*this, s, t, itype);
}

template<typename T>
void ParticleDataImpl<T>::add(const ParticleDataImpl<T>& a)
{
	knPdataAdd<T,T> op( *this, a );
}
template<typename T>
void ParticleDataImpl<T>::sub(const ParticleDataImpl<T>& a)
{
	knPdataSub<T,T> op( *this, a );
}

template<typename T>
void ParticleDataImpl<T>::addConst(const T &s)
{
	knPdataAddScalar<T,T> op( *this, s );
}

template<typename T>
void ParticleDataImpl<T>::addScaled(const ParticleDataImpl<T>& a, const T& factor)
{
	knPdataScaledAdd<T,T> op( *this, a, factor );
}

template<typename T>
void ParticleDataImpl<T>::mult(const ParticleDataImpl<T>& a)
{
	knPdataMult<T,T> op( *this, a );
}

template<typename T>
void ParticleDataImpl<T>::safeDiv(const ParticleDataImpl<T>& a)
{
	knPdataSafeDiv<T> op( *this, a );
}

template<typename T>
void ParticleDataImpl<T>::multConst(const T &s)
{
	knPdataMultScalar<T,T> op( *this, s );
}


template<typename T>
void ParticleDataImpl<T>::clamp(const Real vmin, const Real vmax)
{
	knPdataClamp<T> op( *this, vmin, vmax );
}

template<typename T>
void ParticleDataImpl<T>::clampMin(const Real vmin)
{
	knPdataClampMin<T> op( *this, vmin );
}
template<typename T>
void ParticleDataImpl<T>::clampMax(const Real vmax)
{
	knPdataClampMax<T> op( *this, vmax );
}

template<>
void ParticleDataImpl<Vec3>::clampMin(const Real vmin)
{
	knPdataClampMinVec3 op( *this, vmin );
}
template<>
void ParticleDataImpl<Vec3>::clampMax(const Real vmax)
{
	knPdataClampMaxVec3 op( *this, vmax );
}

template<typename T> KERNEL(pts, reduce=+) returns(T result=T(0.)) T    KnPtsSum(const ParticleDataImpl<T>& val, const ParticleDataImpl<int> *t, const int itype) { if(t && !((*t)[idx]&itype)) return; result += val[idx]; }
template<typename T> KERNEL(pts, reduce=+) returns(Real result=0.) Real KnPtsSumSquare(const ParticleDataImpl<T>& val)    { result += normSquare(val[idx]); }
template<typename T> KERNEL(pts, reduce=+) returns(Real result=0.) Real KnPtsSumMagnitude(const ParticleDataImpl<T>& val) { result += norm(val[idx]); }

template<typename T>
T ParticleDataImpl<T>::sum(const ParticleDataImpl<int> *t, const int itype) const
{
	return KnPtsSum<T>(*this, t, itype);
}
template<typename T>
Real ParticleDataImpl<T>::sumSquare() const
{
	return KnPtsSumSquare<T>(*this);
}
template<typename T>
Real ParticleDataImpl<T>::sumMagnitude() const
{
	return KnPtsSumMagnitude<T>(*this);
}

template<typename T>
KERNEL(pts, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompPdata_Min(const ParticleDataImpl<T>& val)
{
	if(val[idx] < minVal) minVal = val[idx];
}

template<typename T>
KERNEL(pts, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real CompPdata_Max(const ParticleDataImpl<T>& val)
{
	if(val[idx] > maxVal) maxVal = val[idx];
}

template<typename T>
Real ParticleDataImpl<T>::getMin() const
{
	return CompPdata_Min<T>(*this);
}

template<typename T>
Real ParticleDataImpl<T>::getMaxAbs() const
{
	Real amin = CompPdata_Min<T>(*this);
	Real amax = CompPdata_Max<T>(*this);
	return max(fabs(amin), fabs(amax));
}

template<typename T>
Real ParticleDataImpl<T>::getMax() const
{
	return CompPdata_Max<T>(*this);
}

template<typename T>
void ParticleDataImpl<T>::printPdata(IndexInt start, IndexInt stop, bool printIndex)
{
	std::ostringstream sstr;
	IndexInt s = (start>0 ? start : 0                      );
	IndexInt e = (stop>0  ? stop  : (IndexInt)mData.size() );
	s = Manta::clamp(s, (IndexInt)0, (IndexInt)mData.size());
	e = Manta::clamp(e, (IndexInt)0, (IndexInt)mData.size());

	for(IndexInt i=s; i<e; ++i) {
		if(printIndex) sstr << i<<": ";
		sstr<<mData[i]<<" "<<"\n";
	}
	debMsg( sstr.str() , 1 );
}
template<class T> std::string ParticleDataImpl<T>::getDataPointer() {
	std::ostringstream out;
	out << &mData;
	return out.str();
}

// specials for vec3
// work on length values, ie, always positive (in contrast to scalar versions above)

KERNEL(pts, reduce=min) returns(Real minVal=std::numeric_limits<Real>::max())
Real CompPdata_MinVec3(const ParticleDataImpl<Vec3>& val)
{
	const Real s = normSquare(val[idx]);
	if(s < minVal) minVal = s;
}

KERNEL(pts, reduce=max) returns(Real maxVal=-std::numeric_limits<Real>::max())
Real CompPdata_MaxVec3(const ParticleDataImpl<Vec3>& val)
{
	const Real s = normSquare(val[idx]);
	if(s > maxVal) maxVal = s;
}

template<>
Real ParticleDataImpl<Vec3>::getMin() const
{
	return sqrt(CompPdata_MinVec3(*this));
}

template<>
Real ParticleDataImpl<Vec3>::getMaxAbs() const
{
	return sqrt(CompPdata_MaxVec3(*this));  // no minimum necessary here
}

template<>
Real ParticleDataImpl<Vec3>::getMax() const
{
	return sqrt(CompPdata_MaxVec3(*this));
}


// explicit instantiation
template class ParticleDataImpl<int>;
template class ParticleDataImpl<Real>;
template class ParticleDataImpl<Vec3>;

//******************************************************************************
// Implementation
//******************************************************************************

const int DELETE_PART = 20; // chunk size for compression

void ParticleBase::addBuffered(const Vec3& pos, int flag) {
	mNewBufferPos.push_back(pos);
	mNewBufferFlag.push_back(flag);
}
   
template<class S>
void ParticleSystem<S>::clear() {
	mDeleteChunk = mDeletes = 0;
	this->resizeAll(0); // instead of mData.clear
}

template<class S>
IndexInt ParticleSystem<S>::add(const S& data) {
	mData.push_back(data); 
	mDeleteChunk = mData.size() / DELETE_PART;
	this->addAllPdata();
	return mData.size()-1;
}

template<class S>
inline void ParticleSystem<S>::kill(IndexInt idx) {
	assertMsg(idx>=0 && idx<size(), "Index out of bounds");
	mData[idx].flag |= PDELETE; 
	if ( (++mDeletes > mDeleteChunk) && (mAllowCompress) ) compress(); 
}

template<class S>
void ParticleSystem<S>::getPosPdata(ParticleDataImpl<Vec3>& target) const {
	for(IndexInt i=0; i<(IndexInt)this->size(); ++i) {
		target[i] = this->getPos(i);
	}
}
template<class S>
void ParticleSystem<S>::setPosPdata(const ParticleDataImpl<Vec3>& source) {
	for(IndexInt i=0; i<(IndexInt)this->size(); ++i) {
		this->setPos(i, source[i]);
	}
}

template<class S>
void ParticleSystem<S>::transformPositions( Vec3i dimOld, Vec3i dimNew )
{
	const Vec3 factor = calcGridSizeFactor( dimNew, dimOld );
	for(IndexInt i=0; i<(IndexInt)this->size(); ++i) {
		this->setPos(i, this->getPos(i) * factor );
	}
}

// check for deletion/invalid position, otherwise return velocity
KERNEL(pts) returns(std::vector<Vec3> u) template<class S>
std::vector<Vec3> GridAdvectKernel(
	std::vector<S>& p, const MACGrid& vel, const FlagGrid& flags, const Real dt,
	const bool deleteInObstacle, const bool stopInObstacle, const bool skipNew,
	const ParticleDataImpl<int> *ptype, const int exclude)
{
	if ((p[idx].flag & ParticleBase::PDELETE) || (ptype && ((*ptype)[idx] & exclude)) || (skipNew && (p[idx].flag & ParticleBase::PNEW))) {
		u[idx] = 0.; return;
	}
	// special handling
	if(deleteInObstacle || stopInObstacle) {
		if (!flags.isInBounds(p[idx].pos, 1) || flags.isObstacle(p[idx].pos) ) {
			if(stopInObstacle)
				u[idx] = 0.;
			// for simple tracer particles, its convenient to delete particles right away
			// for other sim types, eg flip, we can try to fix positions later on
			if(deleteInObstacle)
				p[idx].flag |= ParticleBase::PDELETE;
			return;
		}
	}
	u[idx] = vel.getInterpolated(p[idx].pos) * dt;
};

// final check after advection to make sure particles haven't escaped
// (similar to particle advection kernel)
KERNEL(pts) template<class S>
void KnDeleteInObstacle(std::vector<S>& p, const FlagGrid& flags) {
	if (p[idx].flag & ParticleBase::PDELETE) return;
	if (!flags.isInBounds(p[idx].pos,1) || flags.isObstacle(p[idx].pos)) {
		p[idx].flag |= ParticleBase::PDELETE;
	} 
}

// try to get closer to actual obstacle boundary
static inline Vec3 bisectBacktracePos(const FlagGrid& flags, const Vec3& oldp, const Vec3& newp)
{
	Real s = 0.;
	for(int i=1; i<5; ++i) {
		Real ds = 1./(Real)(1<<i);
		if (!flags.isObstacle( oldp*(1.-(s+ds)) + newp*(s+ds) )) {
			s += ds;
		}
	}
	return( oldp*(1.-(s)) + newp*(s) );
}

// at least make sure all particles are inside domain
KERNEL(pts) template<class S>
void KnClampPositions(
	std::vector<S>& p, const FlagGrid& flags, ParticleDataImpl<Vec3> *posOld=NULL, bool stopInObstacle=true,
	const ParticleDataImpl<int> *ptype=NULL, const int exclude=0)
{
	if (p[idx].flag & ParticleBase::PDELETE) return;
	if (ptype && ((*ptype)[idx] & exclude)) {
		if(posOld) p[idx].pos = (*posOld)[idx];
		return;
	}
	if (!flags.isInBounds(p[idx].pos,0) ) {
		p[idx].pos = clamp( p[idx].pos, Vec3(0.), toVec3(flags.getSize())-Vec3(1.) );
	} 
	if (stopInObstacle && (flags.isObstacle(p[idx].pos)) ) {
		p[idx].pos = bisectBacktracePos(flags, (*posOld)[idx], p[idx].pos);
	}
}

// advection plugin
template<class S>
void ParticleSystem<S>::advectInGrid(
	const FlagGrid &flags, const MACGrid &vel, const int integrationMode,
	const bool deleteInObstacle, const bool stopInObstacle, const bool skipNew,
	const ParticleDataImpl<int> *ptype, const int exclude)
{
	// position clamp requires old positions, backup
	ParticleDataImpl<Vec3> *posOld = NULL;
	if(!deleteInObstacle) {
		posOld = new ParticleDataImpl<Vec3>(this->getParent());
		posOld->resize(mData.size());
		for(IndexInt i=0; i<(IndexInt)mData.size();++i) (*posOld)[i] = mData[i].pos;
	}

	// update positions
	GridAdvectKernel<S> kernel(mData, vel, flags, getParent()->getDt(), deleteInObstacle, stopInObstacle, skipNew, ptype, exclude);
	integratePointSet(kernel, integrationMode);

	if(!deleteInObstacle) {
		KnClampPositions<S>  (mData, flags, posOld, stopInObstacle, ptype, exclude);
		delete posOld;
	} else {
		KnDeleteInObstacle<S>(mData, flags);
	}
}

KERNEL(pts, single) // no thread-safe random gen yet
template<class S>
void KnProjectParticles(ParticleSystem<S>& part, Grid<Vec3>& gradient) {
	static RandomStream rand (3123984);
	const double jlen = 0.1;
	
	if (part.isActive(idx)) {
		// project along levelset gradient
		Vec3 p = part[idx].pos;
		if (gradient.isInBounds(p)) {
			Vec3 n = gradient.getInterpolated(p);
			Real dist = normalize(n);
			Vec3 dx = n * (-dist + jlen * (1 + rand.getReal()));
			p += dx;            
		}
		// clamp to outer boundaries (+jitter)
		const double jlen = 0.1;
		Vec3 jitter = jlen * rand.getVec3();
		part[idx].pos = clamp(p, Vec3(1,1,1)+jitter, toVec3(gradient.getSize()-1)-jitter);
	}
}

template<class S>
void ParticleSystem<S>::projectOutside(Grid<Vec3> &gradient) {
	KnProjectParticles<S>(*this, gradient);
}

KERNEL(pts) template<class S>
void KnProjectOutOfBnd(ParticleSystem<S> &part, const FlagGrid &flags, const Real bnd, const bool *axis, const ParticleDataImpl<int> *ptype, const int exclude) {
	if(!part.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) return;
	if(axis[0]) part[idx].pos.x = std::max(part[idx].pos.x, bnd);
	if(axis[1]) part[idx].pos.x = std::min(part[idx].pos.x, static_cast<Real>(flags.getSizeX())-bnd);
	if(axis[2]) part[idx].pos.y = std::max(part[idx].pos.y, bnd);
	if(axis[3]) part[idx].pos.y = std::min(part[idx].pos.y, static_cast<Real>(flags.getSizeY())-bnd);
	if(flags.is3D()) {
		if(axis[4]) part[idx].pos.z = std::max(part[idx].pos.z, bnd);
		if(axis[5]) part[idx].pos.z = std::min(part[idx].pos.z, static_cast<Real>(flags.getSizeZ())-bnd);
	}
}

template<class S>
void ParticleSystem<S>::projectOutOfBnd(const FlagGrid &flags, const Real bnd, const std::string &plane, const ParticleDataImpl<int> *ptype, const int exclude) {
	bool axis[6] = { false };
	for(std::string::const_iterator it=plane.begin(); it!=plane.end(); ++it) {
		if(*it=='x') axis[0] = true;
		if(*it=='X') axis[1] = true;
		if(*it=='y') axis[2] = true;
		if(*it=='Y') axis[3] = true;
		if(*it=='z') axis[4] = true;
		if(*it=='Z') axis[5] = true;
	}
	KnProjectOutOfBnd<S>(*this, flags, bnd, axis, ptype, exclude);
}

template<class S>
void ParticleSystem<S>::resizeAll(IndexInt size) {
	// resize all buffers to target size in 1 go
	mData.resize(size);
	for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i)
		mPartData[i]->resize(size);
}

template<class S>
void ParticleSystem<S>::compress() {
	IndexInt nextRead = mData.size();
	for (IndexInt i=0; i<(IndexInt)mData.size(); i++) {
		while ((mData[i].flag & PDELETE) != 0) {
			nextRead--;
			mData[i] = mData[nextRead];
			// ugly, but prevent virtual function calls here:
			for(IndexInt pd=0; pd<(IndexInt)mPdataReal.size(); ++pd) mPdataReal[pd]->copyValue(nextRead, i);
			for(IndexInt pd=0; pd<(IndexInt)mPdataVec3.size(); ++pd) mPdataVec3[pd]->copyValue(nextRead, i);
			for(IndexInt pd=0; pd<(IndexInt)mPdataInt .size(); ++pd) mPdataInt [pd]->copyValue(nextRead, i);
			mData[nextRead].flag = PINVALID;
		}
	}
	if(nextRead<(IndexInt)mData.size()) debMsg("Deleted "<<((IndexInt)mData.size() - nextRead)<<" particles", 1); // debug info

	resizeAll(nextRead);
	mDeletes = 0;
	mDeleteChunk = mData.size() / DELETE_PART;
}

//! insert buffered positions as new particles, update additional particle data
template<class S>
void ParticleSystem<S>::insertBufferedParticles() {
	// clear new flag everywhere
	for(IndexInt i=0; i<(IndexInt)mData.size(); ++i) mData[i].flag &= ~PNEW;

	if(mNewBufferPos.size()==0) return;
	IndexInt newCnt = mData.size();
	resizeAll(newCnt + mNewBufferPos.size());

	for(IndexInt i=0; i<(IndexInt)mNewBufferPos.size(); ++i) {
		int flag = (mNewBufferFlag.size() > 0) ? mNewBufferFlag[i] : 0;
		// note, other fields are not initialized here...
		mData[newCnt].pos  = mNewBufferPos[i];
		mData[newCnt].flag = PNEW | flag;
		// now init pdata fields from associated grids...
		for(IndexInt pd=0; pd<(IndexInt)mPdataReal.size(); ++pd) 
			mPdataReal[pd]->initNewValue(newCnt, mNewBufferPos[i] );
		for(IndexInt pd=0; pd<(IndexInt)mPdataVec3.size(); ++pd) 
			mPdataVec3[pd]->initNewValue(newCnt, mNewBufferPos[i] );
		for(IndexInt pd=0; pd<(IndexInt)mPdataInt.size(); ++pd) 
			mPdataInt[pd]->initNewValue(newCnt, mNewBufferPos[i] );
		newCnt++;
	}
	if(mNewBufferPos.size()>0) debMsg("Added & initialized "<<(IndexInt)mNewBufferPos.size()<<" particles", 2); // debug info
	mNewBufferPos.clear();
	mNewBufferFlag.clear();
}


template<class DATA, class CON>
void ConnectedParticleSystem<DATA,CON>::compress() {
	const IndexInt sz = ParticleSystem<DATA>::size();
	IndexInt *renumber_back = new IndexInt[sz];
	IndexInt *renumber = new IndexInt[sz];
	for (IndexInt i=0; i<sz; i++)
		renumber[i] = renumber_back[i] = -1;
		
	// reorder elements
	std::vector<DATA>& data = ParticleSystem<DATA>::mData;
	IndexInt nextRead = sz;
	for (IndexInt i=0; i<nextRead; i++) {
		if ((data[i].flag & ParticleBase::PDELETE) != 0) {
			nextRead--;
			data[i] = data[nextRead];
			data[nextRead].flag = 0;           
			renumber_back[i] = nextRead;
		} else 
			renumber_back[i] = i;
	}
	
	// acceleration structure
	for (IndexInt i=0; i<nextRead; i++)
		renumber[renumber_back[i]] = i;
	
	// rename indices in filaments
	for (IndexInt i=0; i<(IndexInt)mSegments.size(); i++)
		mSegments[i].renumber(renumber);
		
	ParticleSystem<DATA>::mData.resize(nextRead);
	ParticleSystem<DATA>::mDeletes = 0;
	ParticleSystem<DATA>::mDeleteChunk = ParticleSystem<DATA>::size() / DELETE_PART;
	
	delete[] renumber;
	delete[] renumber_back;
}

template<class S>
ParticleBase* ParticleSystem<S>::clone() {
	ParticleSystem<S>* nm = new ParticleSystem<S>(getParent());
	if(this->mAllowCompress) compress();
	
	nm->mData = mData;
	nm->setName(getName());
	this->cloneParticleData(nm);
	return nm;
}

template<class DATA,class CON>
ParticleBase* ConnectedParticleSystem<DATA,CON>::clone() {
	ConnectedParticleSystem<DATA,CON>* nm = new ConnectedParticleSystem<DATA,CON>(this->getParent());
	if(this->mAllowCompress) compress();
	
	nm->mData = this->mData;
	nm->mSegments = mSegments;
	nm->setName(this->getName());
	this->cloneParticleData(nm);
	return nm;
}

template<class S>  
std::string ParticleSystem<S>::infoString() const { 
	std::stringstream s;
	s << "ParticleSys '" << getName() << "'\n-> ";
	if(this->getNumPdata()>0) s<< "pdata: "<< this->getNumPdata();
	s << "parts: " << size();
	//for(IndexInt i=0; i<(IndexInt)mPartData.size(); ++i) { sstr << i<<":" << mPartData[i]->size() <<" "; }
	return s.str();
}
	
template<class S>  
inline void ParticleSystem<S>::checkPartIndex(IndexInt idx) const {
	IndexInt mySize = this->size();
	if (idx<0 || idx > mySize ) {
		errMsg( "ParticleBase " << " size " << mySize << " : index " << idx << " out of bound " );
	}
}
	
inline void ParticleDataBase::checkPartIndex(IndexInt idx) const {
	IndexInt mySize = this->getSizeSlow();
	if (idx<0 || idx > mySize ) {
		errMsg( "ParticleData " << " size " << mySize << " : index " << idx << " out of bound " );
	}
	if ( mpParticleSys && mpParticleSys->getSizeSlow()!=mySize ) {
		errMsg( "ParticleData " << " size " << mySize << " does not match parent! (" << mpParticleSys->getSizeSlow() << ") " );
	}
}

// set contents to zero, as for a grid
template<class T>
void ParticleDataImpl<T>::clear() {
	for(IndexInt i=0; i<(IndexInt)mData.size(); ++i) mData[i] = 0.;
}

template class ParticleSystem<BasicParticleData>;
template class ParticleSystem<TurbulenceParticleData>;
template class ParticleSystem<VortexParticleData>;
} // namespace
