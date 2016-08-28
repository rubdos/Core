Summary:            A free montecarlo raytracing engine with plugins for Blender and Maya3D
Name:               thebounty
Version:            0.1.6.rc4
Release:            1%{?dist}
License:            LGPLv3+
Group:              Applications/Multimedia
Source0: 			https://github.com/TheBounty/Core/archive/%{version}.tar.gz
URL:                https://www.thebountyrenderer.org/
ExcludeArch:        s390 s390x
BuildRequires:      cmake swig
BuildRequires:      zlib-devel,libpng-devel,OpenEXR-libs,OpenEXR-devel
BuildRequires:      ilmbase-devel freetype-devel libxml2-devel libjpeg-devel
BuildRequires:      libtiff-devel python%{python3_pkgversion}-devel libstdc++-devel
Requires:           zlib,libpng,OpenEXR-libs,OpenEXR
Requires:           ilmbase freetype libxml2 libjpeg
Requires:           libtiff, python(abi)>=3.4, libstdc++, system-python-libs

%package devel
Summary: 			TheBounty development files and headers
Group:              Applications/Multimedia
Provides: 			thebounty-devel

%description
TheBounty is a free montecarlo raytracing engine.

Raytracing is a rendering technique for generating realistic images by tracing the path of light through a 3D scene. A render engine consists of a specialised computer program that interacts with a host 3D application to provide specific raytracing capabilities "on demand".

The TheBounty engine is supported in Blender and Wings 3D.

%description devel
TheBounty is a free montecarlo raytracing engine. This package contains the headers for building TheBounty plugins.

Raytracing is a rendering technique for generating realistic images by tracing the path of light through a 3D scene. A render engine consists of a specialised computer program that interacts with a host 3D application to provide specific raytracing capabilities "on demand".

The TheBounty engine is supported in Blender and Wings 3D.

%prep
ls
%autosetup -n %{name}-%{version}

%build
mkdir build
pushd build
%define pythonversion %(python3 --version | awk '{print $2}' | awk -F . '{print $1"."$2}')
cmake -DWITH_QT=off -DYAF_PY_VERSION=%{pythonversion} \
	-DWITH_YAF_PY_BINDINGS=ON \
	-DYAF_BIN_DIR=%{_bindir} -DYAF_LIB_DIR=%{_libdir} -DYAF_PLUGIN_DIR=%{_libdir}/thebounty/ \
	-DYAF_BINDINGS_PY_DIR=%{_libdir}/python%{pythonversion}/site-packages/ \
	-DBUILDRELEASE=ON \
	-DCMAKE_BUILD_TYPE=Release ..
%make_build
popd

%install
pushd build

%make_install
popd

mkdir -p %{buildroot}/%{_includedir}/thebounty/
install build/*.h %{buildroot}/%{_includedir}/thebounty/
find include -iname \*.h -exec install "{}" %{buildroot}/%{_includedir}/thebounty/ \;


%files
%doc README.md
%{_bindir}/*
%{_libdir}/*

%files devel
%{_bindir}/*
%{_libdir}/*
%{_includedir}/*

%changelog
* Sat Aug 27 2016 Ruben De Smet <ruben.de.smet@thebountyrenderer.org> 0.1.6.rc4-1
- First automatic build
