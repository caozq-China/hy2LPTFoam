#include "dsmcCloud.H"
#include "CloudFunctionObject.H"
//#include "FacePostProcessing.H"
#include "dsmcFaceTrackerFO.H"

#include "runTimeSelectionTables.H"
/*
    defineNamedTemplateTypeNameAndDebug
    (
        Foam::CloudFunctionObject<Foam::dsmcCloud>,
        0
    );
    namespace Foam
    {
        defineTemplateRunTimeSelectionTable
        (
            CloudFunctionObject<dsmcCloud>,
            dictionary
        );
    }
*/
makeCloudFunctionObject(dsmcCloud);

namespace Foam
{
//makeCloudFunctionObjectType(FacePostProcessing, dsmcCloud);
//makeCloudFunctionObjectTypeNT(dsmcFaceTrackerFO, dsmcCloud);
}

